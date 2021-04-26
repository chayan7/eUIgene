
Steps:

1. Download and Install Conda (https://docs.anaconda.com/anaconda/install/)
2. Download the eUIgene tool from github
3. Install Conda Environment:
	cd eUIgene/
	bash ./build.sh
4. Activate eUIgene Environment
5. Run eUIgene script: # eUIgene.py is located inside script directory
	python3 eUIgene.py -ot exampleInput_flagsout_TreeOrder_operon.tsv -p 300 -o exampleOutput
