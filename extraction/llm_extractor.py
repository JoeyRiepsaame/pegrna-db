"""LLM-based extraction using Claude or OpenAI-compatible APIs."""
import json
import re
from typing import Optional

from rich.console import Console

import config
from extraction.schemas import PegRNAExtracted

console = Console()

EXTRACTION_PROMPT = """You are an expert in CRISPR prime editing. Extract all experimentally validated pegRNA and epegRNA entries from the following scientific paper text/table data.

For EACH pegRNA/epegRNA entry, extract these fields (use null if not found):
- entry_name: identifier used in the paper (e.g., "pegRNA-1", "HEK3-pegRNA")
- pegrna_type: "pegRNA" or "epegRNA"
- spacer_sequence: the ~20nt guide/spacer sequence (DNA, no scaffold)
- pbs_sequence: primer binding site sequence
- pbs_length: PBS length in nucleotides
- rtt_sequence: reverse transcriptase template sequence
- rtt_length: RTT length in nucleotides
- three_prime_extension: 3' structural motif name (e.g., "evopreQ1", "mpknot", "tevopreQ1")
- linker_sequence: linker between extension and 3' motif
- full_sequence: complete pegRNA sequence if available
- nicking_sgrna_seq: nicking sgRNA sequence for PE3/PE3b
- target_gene: gene name (e.g., "HEK3", "FANCF", "EMX1")
- target_locus: genomic coordinates if available
- target_organism: species (e.g., "human", "mouse")
- edit_type: "substitution", "insertion", or "deletion"
- edit_description: description of the edit (e.g., "C>T at +5")
- intended_mutation: specific nucleotide change
- prime_editor: which PE system (e.g., "PE2", "PE3", "PE3b", "PE4max", "PE5max", "PE7")
- cell_type: cell line or type (e.g., "HEK293T", "HeLa", "K562")
- delivery_method: how pegRNA was delivered (e.g., "plasmid", "RNP", "mRNA")
- editing_efficiency: percentage editing efficiency (number only, 0-100)
- product_purity: percentage of correct edits among all edits (0-100)
- indel_frequency: percentage of indels (0-100)

IMPORTANT:
- Only include entries with EXPERIMENTAL data (actual editing efficiency measurements)
- Include ALL entries from tables, even if there are hundreds
- Sequences should be DNA (ACGT), convert U to T
- Efficiency values should be numbers (not strings), as percentages 0-100
- If a table has multiple conditions for the same pegRNA, create separate entries

Return a JSON array of objects. Nothing else."""


def _call_anthropic(prompt: str) -> Optional[str]:
    """Call Claude API."""
    try:
        from anthropic import Anthropic
        client = Anthropic(api_key=config.ANTHROPIC_API_KEY)
        message = client.messages.create(
            model=config.LLM_MODEL,
            max_tokens=config.LLM_MAX_TOKENS,
            messages=[{"role": "user", "content": prompt}],
        )
        console.print(
            f"(Claude: {message.usage.input_tokens}in / {message.usage.output_tokens}out)"
        )
        return message.content[0].text.strip()
    except Exception as e:
        console.print(f"[yellow]Claude API failed: {e}[/yellow]")
        return None


def _call_openai_compatible(prompt: str) -> Optional[str]:
    """Call OpenAI-compatible API (OpenRouter, DeepSeek, etc.)."""
    import requests

    # Try OpenRouter first (has many models)
    openrouter_key = config.OPENROUTER_API_KEY
    if openrouter_key:
        try:
            resp = requests.post(
                "https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {openrouter_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "anthropic/claude-sonnet-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "max_tokens": config.LLM_MAX_TOKENS,
                },
                timeout=120,
            )
            if resp.status_code == 200:
                data = resp.json()
                return data["choices"][0]["message"]["content"].strip()
            else:
                console.print(f"[yellow]OpenRouter failed: {resp.status_code}[/yellow]")
        except Exception as e:
            console.print(f"[yellow]OpenRouter error: {e}[/yellow]")

    # Try ChatGPT
    chatgpt_key = config.CHATGPT_API_KEY
    if chatgpt_key:
        try:
            resp = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {chatgpt_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "gpt-4o",
                    "messages": [{"role": "user", "content": prompt}],
                    "max_tokens": config.LLM_MAX_TOKENS,
                },
                timeout=120,
            )
            if resp.status_code == 200:
                data = resp.json()
                return data["choices"][0]["message"]["content"].strip()
            else:
                console.print(f"[yellow]ChatGPT failed: {resp.status_code}[/yellow]")
        except Exception as e:
            console.print(f"[yellow]ChatGPT error: {e}[/yellow]")

    return None


def extract_with_llm(
    paper_text: str,
    paper_title: str = "",
    paper_abstract: str = "",
    max_text_length: int = 50000,
) -> list[PegRNAExtracted]:
    """Extract pegRNA data from paper text using LLM API.

    Tries Claude first, then falls back to OpenRouter/ChatGPT.
    """
    # Truncate if needed
    if len(paper_text) > max_text_length:
        paper_text = paper_text[:max_text_length] + "\n[...truncated...]"

    context = f"Paper title: {paper_title}\n"
    if paper_abstract:
        context += f"Abstract: {paper_abstract}\n"
    context += f"\n---\n\n{paper_text}"

    full_prompt = f"{EXTRACTION_PROMPT}\n\n---\n\n{context}"

    console.print("[bold]Sending to LLM for extraction...[/bold]")

    # Try Claude first
    response_text = None
    if config.ANTHROPIC_API_KEY:
        response_text = _call_anthropic(full_prompt)

    # Fallback to OpenAI-compatible APIs
    if not response_text:
        response_text = _call_openai_compatible(full_prompt)

    if not response_text:
        console.print("[red]All LLM backends failed[/red]")
        return []

    entries = _parse_llm_response(response_text)
    console.print(f"[green]LLM extracted {len(entries)} entries[/green]")
    return entries


def extract_table_with_llm(
    table_text: str,
    paper_title: str = "",
) -> list[PegRNAExtracted]:
    """Extract from a specific table that rule-based extraction couldn't handle."""
    prompt = (
        "Extract all pegRNA/epegRNA entries from this table. "
        "The table may have non-standard column names or formatting.\n\n"
        f"Paper: {paper_title}\n\n"
        f"Table data:\n{table_text}"
    )
    return extract_with_llm(prompt, paper_title)


def _parse_llm_response(response: str) -> list[PegRNAExtracted]:
    """Parse LLM JSON response into validated PegRNAExtracted objects."""
    json_str = response

    # Handle markdown code blocks
    if "```" in response:
        match = re.search(r"```(?:json)?\s*\n?(.*?)```", response, re.DOTALL)
        if match:
            json_str = match.group(1).strip()

    try:
        data = json.loads(json_str)
    except json.JSONDecodeError:
        start = response.find("[")
        end = response.rfind("]")
        if start >= 0 and end > start:
            try:
                data = json.loads(response[start : end + 1])
            except json.JSONDecodeError:
                console.print("[red]Failed to parse LLM response as JSON[/red]")
                return []
        else:
            console.print("[red]No JSON array found in LLM response[/red]")
            return []

    if not isinstance(data, list):
        data = [data]

    entries = []
    for item in data:
        try:
            item["confidence_score"] = 0.6
            entry = PegRNAExtracted(**item)
            if entry.spacer_sequence or entry.target_gene or entry.full_sequence:
                entries.append(entry)
        except Exception as e:
            console.print(f"[yellow]Skipping invalid LLM entry: {e}[/yellow]")

    return entries
