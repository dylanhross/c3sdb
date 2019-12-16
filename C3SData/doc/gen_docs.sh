#!/bin/bash

python3 -m pydoc -w C3SData
python3 -m pydoc -w C3SData.data
python3 -m pydoc -w C3SData.db_interface

mv -f C3SData*.html C3SData/doc/
