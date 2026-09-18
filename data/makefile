all: 
	-cd swap && ./swap420
	cd swap/output_swap && python3 pressure_bc.py
	cd swap/output_swap && python3 pressure_ini.py
	cd swap/output_swap && python3 water_content.py

input:
	cd swap && python3 input_geology.py

clean:
	rm -f ../src/bottom_pressure
	rm -f ../src/pressure_ini
	rm -f ../src/top_pressure
	rm -f ../src/wc_ini
	rm -f ._*

