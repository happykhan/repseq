"""Rich console and messaging helpers (bactscout-style)."""

from __future__ import annotations

from rich.console import Console

console = Console()

_STYLES: dict[str, dict[str, str]] = {
    "error":   {"style": "bold red",    "emoji": "❌",  "prefix": "ERROR"},
    "warning": {"style": "bold yellow", "emoji": "⚠️ ", "prefix": "WARNING"},
    "success": {"style": "bold green",  "emoji": "✅",  "prefix": "SUCCESS"},
    "info":    {"style": "bold blue",   "emoji": "ℹ️ ", "prefix": "INFO"},
    "debug":   {"style": "dim cyan",    "emoji": "🔍",  "prefix": "DEBUG"},
}


def print_message(message: str, msg_type: str = "info") -> None:
    """Print a coloured, prefixed console message using Rich."""
    cfg = _STYLES.get(msg_type.lower(), _STYLES["info"])
    emoji  = f"{cfg['emoji']} "
    prefix = f"[{cfg['style']}]{cfg['prefix']}[/]"
    body   = f"[{cfg['style']}]{message}[/]"
    console.print(f"{emoji}{prefix}: {body}")


def print_header(title: str) -> None:
    """Print a bold magenta section header."""
    console.print(f"\n[bold magenta]{'=' * 60}[/]")
    console.print(f"[bold magenta]{title.center(60)}[/]")
    console.print(f"[bold magenta]{'=' * 60}[/]\n")
