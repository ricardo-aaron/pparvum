#!/usr/bin/perl -w

# $dir is the path to the fastq files
$dir = "/home/rchavez/TexasTech/Mary/fastq/";

# $diro is the path for the output directory
$diro = "/castle/rchavez/Mary/kallistos/kallisto_Mary_fastqs_filtered/";
system("mkdir ".$diro) if (! -d $diro);

# index is the kallisto index file, with its path
$index = "/castle/rchavez/Mary/transcriptomes/transcriptome__drap_ALL__Mary_fastqs/final_tpm1_Mary/index_kallisto/Ppar_transcripts";
$cores = "8";

%hash = ();
$dh = undef;
opendir($dh, $dir) or die $!;
while ($file = readdir($dh))
{
	undef $sample;
	if ($file =~ /(\w)_S\d{1,2}_L001_R[12]_\d{3}\.fastq/)
	{
		$sample = $1;
		$hash{$sample}{$file} = 1;
	}
}
closedir($dh);

foreach $sample (sort(keys(%hash)))
{
	$diro_sample = $diro.$sample."/";
	#next if (-f $diro_sample."abundance.tsv");
	#
	undef $sh;
	$sh = "kallisto quant".
		" -i ".$index.
		" -o ".$diro_sample.
		#" --single".
		#" -l 300".
		#" -s 80" .
		" --bias".
		#" -b 100".
		" --seed 777".
		" -t ".$cores.
		" ".$dir.$a[0].
		" ".$dir.$a[1].
		" 2> ".$diro.$sample.".log";
	print $sh."\n\n";
	system($sh);
}
