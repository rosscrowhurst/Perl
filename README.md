# Perl

Miscellaneous perl code for Bioinformatics

## Genewise Prediction

Simple BASH workflow to take TBLASTn comparison results, set up up genewsie commands, run then parse and finally run evigene

```{bash}
# Setup Blast Database and run blast comparison
PROJECTNAME=CompWithRosidProts
CPUS=16
WORKING_DIR=/workspace
REFERENCE_DIR=/reference
QUERY_DIR=/TBLASTNRosidProteins
REFERENCE_NAME=reference
QUERY_NAME=Rosid_Predicted_Proteins_unique
QUERY_FASTA=${QUERY_DIR}/${QUERY_NAME}.faa
REFERENCE_FASTA=${REFERENCE_DIR}/${REFERENCE_NAME}.fa
BLAST_DB_DIR=${WORKING_DIR}/Reference/BLAST_IDX
mkdir -p ${BLAST_DB_DIR}
cd ${BLAST_DB_DIR}
ln -s ${REFERENCE_FASTA} ${REFERENCE_NAME}
formatdb -i ${REFERENCE_NAME} -p F -o T
echo "blastall -p tblastn -m 8 -a $CPUS -i ${QUERY_FASTA} -d ${BLAST_DB_DIR}/${REFERENCE_NAME} -o ${QUERY_NAME}_${REFERENCE_NAME}_genes.tblastn" > run_tblastn.sh
bsub -n $CPUS -m "${PROJECTNAME}" -o `pwd`/ol.tblastn.stdout -e `pwd`/ol.tblastn.stderr < run_tblastn.sh

# Explode reference and query fasta file to get individual files (needed for genewise commands later)
REF_INDIVID_FILES_DIR=${REFERENCE_DIR}/Genome/Fasta/IndividualScaffolds
mkdir -p ${REF_INDIVID_FILES_DIR}
cp ${REFERENCE_FASTA} ${REF_INDIVID_FILES_DIR}
cd ${REF_INDIVID_FILES_DIR}
fastaexplode ${REFERENCE_FASTA}
rm ${REFERENCE_FASTA}
QUERY_INDIVID_FILES_DIR=QUERY_DIR/Protein/Fasta/IndividuaProteins
mkdir -p ${QUERY_INDIVID_FILES_DIR}
cp ${QUERY_FASTA} ${QUERY_INDIVID_FILES_DIR}
cd ${QUERY_INDIVID_FILES_DIR}
fastaexplode ${QUERY_FASTA}
rm ${QUERY_FASTA}

# Run tabular blast parser on tblastn results
GENEWISE_DIR=${WORKING_DIR}/GeneWise
blast92gff3_rnc.pl -gff3Outfile=modelproteins-mygenome.gff3 \
 -proteinFastaFile=${QUERY_FASTA} \
 -tblastnFile=${QUERY_NAME}_${REFERENCE_NAME}_genes.tblastn \
 -refFastaDir=${REF_INDIVID_FILES_DIR} \
 -proteinFastaDir=${QUERY_INDIVID_FILES_DIR} \
 -genewiseResultsOutDir=${GENEWISE_DIR}

# Split genewise commands from above parsing - produces files named like myfile_1
cd ${GENEWISE_DIR}
file_splitter_by_line.pl run_genewise.sh 40

# Run genewise commands directly or submit to job schedular (not shown)
for n in `seq 0 39`; do (sh myfile_$n 1>myfile_$n.out 2>myfile_$n.err &); done

# Get list of individual results files from genewise runs
ls | grep ".stdout"  > genewise.out.list
for file in `cat genewise.out.list`; do cat $file >> GeneWise.predictions; done

# Parse the genewise predictions
parse_genewise_prediction_results.pl -g=GeneWise.predictions

# Run evigene
tr2aacds.pl -NCPU 16 -MAXMEM 100000 -logfile `pwd`/evigene.log -cdnaseq `pwd`/GeneWise.predictions.fna -MINCDS 90 1>evigene.stderr 2>evigene.stdout &
