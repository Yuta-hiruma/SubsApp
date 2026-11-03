all: 
	python3 run_benchmark.py

plot:
	python3 plotgmt.py
	bash plot_Se_vG.sh
	
clean:
	rm -f output/gmt/*
	rm -f output/*.dat
	rm -f output/*.png
	rm -f output/._*


very_clean:
	rm -f job.*
	rm -f output/*
	rm -f output/._*
	rm -f output/gmt/*
	rm -f ._*


