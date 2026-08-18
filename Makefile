.PHONY: setup pipeline dashboard test

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python run_pipeline.py

dashboard:
	STREAMLIT_BROWSER_GATHER_USAGE_STATS=false streamlit run dashboard.py --server.headless true

test:
	pytest -q -p no:cacheprovider
