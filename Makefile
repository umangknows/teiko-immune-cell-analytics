.PHONY: setup pipeline dashboard test

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python run_pipeline.py

dashboard:
	STREAMLIT_BROWSER_GATHER_USAGE_STATS=false python -m streamlit run dashboard.py --server.headless true

test:
	python -m pytest -q -p no:cacheprovider
