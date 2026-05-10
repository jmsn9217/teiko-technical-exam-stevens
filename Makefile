PYTHON ?= python

.PHONY: setup pipeline dashboard clean

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) -m src.analysis

dashboard:
	$(PYTHON) -m streamlit run dashboard/app.py

clean:
	rm -f teiko.db
	rm -rf outputs
