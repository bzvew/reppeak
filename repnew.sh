#! /bin/bash
SECONDS=0
DIR=$( cd $(dirname $0) ; pwd -P )
printf "script path = %s\n" $DIR
#source ${DIR}/generateFasta.sh
# set flog vars to empty
genomeSize= file= temp= cor=false outputfile="repnew.log" canuPath="" faidxPath="" javaPath="" minOverlapLength=500 minReadLength=1000 
while getopts f:o:r:s:e:h:n:l:t:p:j:b:q:w:d:x:c:g:a:m:z:u: opt
do
	case $opt in
		f)	file=$OPTARG
			;;
		s)	genomeSize=$OPTARG
			;;
		h)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))
			continue
			fi
			maxThreads=$OPTARG
			;;
		e)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))
			continue
			fi
			maxMemory=$OPTARG
			;;
		r)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))
			continue
			fi
			minReadLength=$OPTARG
			;;
		o)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))
			continue
			fi
			minOverlapLength=$OPTARG
			;;
		j)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))  
			continue
			fi
			javaPath=$OPTARG
			;;
		t)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))  
			continue
			fi
			temp=$OPTARG
			;;
		a)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))  
			continue
			fi
			faidxPath=$OPTARG
			;;
		c)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))  
			continue
			fi
			cor=$OPTARG
			;;
		z)	if [[ $OPTARG = -* ]]; then
			((OPTIND--))  
			continue
			fi
			outputfile=$OPTARG
	esac
done
shift $((OPTIND - 1))
printf "file=%s\n" $file
printf "genomeSize=%s\n" $genomeSize
#printf "range=%s\n"	$range
printf "correction=%s\n" $cor
#printf "lines=%s\n" $lines
printf "temp folder=%s\n" ${temp}
#printf "breaks=%d\n" ${breaks}
#printf "n1=%d\n" $n1
#printf "n2=%d\n" $n2
#printf "degree=%d\n" $degree
#printf "commu_size=%d\n" $commu_size
#printf "drops=%d\n" $drops
#printf "ratio=%s\n" $ratio
#printf "window=%d\n" $window
#printf "outputfile=%s\n" $outputfile
home=$(pwd)
orifile="$(cd "$(dirname "$file")" && pwd)/$(basename "$file")"
printf "original place=%s\n" $home
#mkdir $temp
if [ -z $javaPath ]
then
	javaPath=$(command -v java)	
	javaPath=${javaPath%java}
else
	ORGPATH=`pwd`
	RELPATH=$javaPath
	cd $RELPATH
	javaPath=`pwd`
	cd $ORGPATH
fi
#if [ -z $canuPath ]
#then
#	canuPath=$(command -v canu)	
#	canuPath=${canuPath%canu}
#else
#	oRGPATH=`pwd`
#	rELPATH=$canuPath
#	cd $RELPATH
#	canuPath=`pwd`
#	cd $ORGPATH
#fi
canupath="canu-1.4/Linux-amd64/bin"
canupath="$DIR/$canuPath"
if [ -z $faidxPath ]
then
	faidxPath=$(command -v faidx)	
	faidxPath=${faidxPath%faidx}
else
	ORGPATH=`pwd`
	RELPATH=$faidxPath
	cd $RELPATH
	faidxPath=`pwd`
	cd $ORGPATH
fi

printf "canu path is %s\n" $canuPath
printf "min Read Length is %s\n" $minReadLength
printf "min Overlap Length is %s\n" $minOverlapLength
printf "faidx path is %s\n" $faidxPath
printf "java path is %s\n" $javaPath
if [ $cor = true ]
then
	
	printf "Use raw reads\n"
	#$canuPath/canu -p "assem" -d $temp genomeSize="$genomeSize"  saveReadCorrections=T maxThreads=$maxThreads maxMemory=$maxMemory java=$javaPath/java corOutCoverage=400 gnuplotTested=true minReadLength=$minReadLength minOverlapLength=$minOverlapLength corMinCoverage=0 -pacbio-raw "$file"
	printf "the folder is %s\n" $temp
	cd $temp
	printf "process reads\n"
	#bwa index assem.contigs.fasta
#	bwa mem -x pacbio assem.contigs.fasta -a $orifile > reads_assem.sam
#	samtools view -hb  reads_assem.sam | samtools sort - > reads_assem_sorted.bam
#	bedtools genomecov -ibam reads_assem_sorted.bam -bg > reads_assem_sorted.bedgraph
	awk '{print $2,$3,$4>$1}' reads_assem_sorted.bedgraph
	parallel --eta -j 36 --load 80% --noswap 'lines=$(wc -l {} | cut -d " " -f1); python3.4 ~/repnew/extract.py {} $lines 2 30 63' ::: tig*
	cat $(find . -name "tig*.bed" -size +0) > all_tig.bed
	faidx -b all_tig.bed -l assem.contigs.fasta > result.fa
else
	printf "Use corrected reads\n"
	#$canuPath/canu -p "assem" -d $temp genomeSize="$genomeSize" corOutCoverage=400 java=$javaPath/java gnuplotTested=true maxThreads=$maxThreads maxMemory=$maxMemory corMinCoverage=0 minReadLength=$minReadLength minOverlapLength=$minOverlapLength -pacbio-corrected "$file"
	printf "the folder is %s\n" $temp
	cd $temp
	printf "process reads\n"
	#bwa index assem.contigs.fasta
#	bwa mem -x pacbio assem.contigs.fasta -a $orifile > reads_assem.sam
	#samtools view -hb  reads_assem.sam | samtools sort - > reads_assem_sorted.bam
#	bedtools genomecov -ibam reads_assem_sorted.bam -bg > reads_assem_sorted.bedgraph
	awk '{print $2,$3,$4>$1}' reads_assem_sorted.bedgraph
	parallel --eta -j 36 --load 80% --noswap 'lines=$(wc -l {} | cut -d " " -f1); python3.4 ~/repnew/extract.py {} $lines 2 30 63' ::: tig*
	cat $(find . -name "tig*.bed" -size +0) > all_tig.bed
	faidx -b all_tig.bed -l assem.contigs.fasta >result.fa 
fi	
echo $parameters >> ${home}/${outputfile}
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." >> ${home}/${outputfile}
~/replong/dro_evaluation.sh
#rm -rf $temp
cd $home
