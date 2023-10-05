python ../regiosqm/regiosqm.py -g ./example.smiles > ./example.csv

foreach ($file in Get-ChildItem -Filter "*mop") { mopac $file }

python ../regiosqm/regiosqm.py -a ./example.smiles ./example.csv > ./example_results.csv