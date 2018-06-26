all:
	pandoc --title-prefix "Kraken 2 Manual" \
	       --include-in-header head.html \
	       --include-before-body top.html \
	       --from markdown --to html \
	       --table-of-contents \
	       --css kraken.css \
	       --output MANUAL.html \
	       < MANUAL.markdown
