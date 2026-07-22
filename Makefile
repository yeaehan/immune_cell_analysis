.PHONY: setup pipeline dashboard

setup:
	python -m pip install -e .

pipeline:
	python load_data.py
	MPLCONFIGDIR=.mplconfig python analysis.py

dashboard:
	streamlit run dashboard.py
