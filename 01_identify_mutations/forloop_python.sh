for file in *_new_8species.fasta; do file="${file%_new_8species.fasta}"; echo $file; python mutation.py $file; done
