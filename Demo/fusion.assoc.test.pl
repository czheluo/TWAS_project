#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Math::Round;
my ($gwas,$ss,$bim,$pos,$plot,$wp,$ldchr,$chrlist,$fout);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"GWAS:s"=>\$gwas,
	"sumstats:s"=>\$ss,
	"refldchr:s"=>\$ldchr,
	#"out:s"=>\$fout,
	"weightdir:s"=>\$wp,
	"pos:s"=>\$pos,
	"bim:s"=>\$bim,
    "plot:s"=>\$plot,
	"chrlist:s"=>\$chrlist,
			) or &USAGE;
&USAGE unless ($gwas);

#mkdir $fout if (!-d $fout);
#mkdir $wp if (!-d $wp);
$wp=ABSOLUTE_DIR($wp);
#$ldchr=ABSOLUTE_DIR($ldchr);
my $assoc="/mnt/ilustre/centos7users/meng.luo/project/RNA/linannan_eQTL/newGENOME/fusion_twas/fusion_twas/FUSION.assoc_test.R";

open GWAS,$gwas;
my %ch;
while (<GWAS>) {
	chomp;
	next if(/^SNP/);
	my @all=split/\s+/,$_;
	my $zscore=round($all[9]/$all[10]*1000)/1000;
	$ch{$all[0]}=$zscore;
}
close GWAS;
open BM,$bim;
open Out,">$ss";
print Out "SNP\tA1\tA1\tZ\n";
while (<BM>) {
	chomp;
	my(undef,$snp,undef,undef,$a1,$a2)=split/\s+/,$_;
	if (exists $ch{$snp}) {
		print Out "$snp\t$a1\t$a2\t$ch{$snp}\n";
	}	
}
close Out;
close BM;
open CS,$chrlist;
while (<CS>) {
	chomp;
	my(undef,$chr)=split/\s+/,$_;
	$chr=~s/chr//g;
    `Rscript $assoc --sumstats $ss  --weights $pos --weights_dir $wp --ref_ld_chr $ldchr/chr$chr --chr $chr --out $chr.dat`;
}
close CS;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:

	eg: perl $Script -GWAS N.FaSTLMM.txt -sumstats N.sumstats -refldchr /mnt/ilustre/centos7users/meng.luo/Pipeline/locuszoom/locuszoom/data/1000G/genotypes/2021-04-08/napus/ -weightdir ./ -pos Bna.pos -bim pop.bim -chrlist chr.list
	
Usage:
  Options:
	"GWAS:s"=>\$gwas,
	"sumstats:s"=>\$ss,
	"refldchr:s"=>\$ldchr,
	"weightdir:s"=>\$wp,
	"pos:s"=>\$pos,
	"bim:s"=>\$bim,
    "plot:s"=>\$plot,
	"chrlist:s"=>\$chrlist,

USAGE
        print $usage;
        exit;
}
