# Agent Instructions

This repository is a Quarto website managed with `uv`.

## Quarto

- Use `uv` to run Quarto commands so the project uses the pinned Python environment and the `quarto-cli` dependency from `pyproject.toml`.
- Do not assume a globally installed `quarto` binary is available.
- Prefer these commands:
  - `uv run quarto render`
  - `uv run quarto preview`
  - `uv run quarto check`
- If dependencies are missing, run `uv sync` before running Quarto.
- Do *not* run these commands unless explicitly asked to, or if you're debugging.

## Project Notes

- The Quarto configuration is in `_quarto.yml`.
- Generated site output lives in `_site/`.
- Preserve existing content and generated files unless the task explicitly requires changing them.
