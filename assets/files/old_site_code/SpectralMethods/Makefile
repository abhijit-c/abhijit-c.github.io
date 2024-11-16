default:
	pandoc -s main.tex --mathjax --highlight-style pygments -o temp.html
	echo -e '---\nlayout: post\ntitle: "Spectral Methods and the Heat Equation"\ndate: 2019-06-12\n---' | cat - temp.html > temp
	mv temp ../../_posts/2019-06-12-Spectral-Methods.html
	rm temp.html
