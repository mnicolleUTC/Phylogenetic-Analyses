#!/bin/bash
echo 'Checking for pyenv version'
pyenv --version
echo 'Loading python 3.9'
pyenv local 3.9
echo 'Installing requirements'
pip install -r requirements.txt
echo 'Run main.py'
python main.py
