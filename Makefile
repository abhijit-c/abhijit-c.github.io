QUARTO := uv run quarto

.PHONY: help preview publish

preview:
	$(QUARTO) preview

help:
	@echo "make preview  - serve the site locally with live reload"
