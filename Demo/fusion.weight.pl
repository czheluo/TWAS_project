#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bfile,$exp,$pos,$flat,$model,$wp,$fout);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"expression:s"=>\$exp,
	"bfile:s"=>\$bfile,
	"out:s"=>\$fout,
	"models:s"=>\$model,
	"weight:s"=>\$wp,
	"pos:s"=>\$pos,
	"refflat:s"=>\$flat,
			) or &USAGE;
&USAGE unless ($fout);
mkdir $fout if (!-d $fout);
mkdir $wp if (!-d $wp);
$wp=ABSOLUTE_DIR($wp);
$fout=ABSOLUTE_DIR($fout);
my $models||="enet";
my $weight="/mnt/ilustre/centos7users/meng.luo/project/RNA/linannan_eQTL/newGENOME/fusion_twas/fusion_twas/examples/FUSION.compute_weights.R";
my $gcta="/mnt/ilustre/centos7users/dna/.env/bin/gcta";
my $gemma="/mnt/ilustre/centos7users/meng.luo/Pipeline/DNA/09.GWAS-SNP/gemma-0.94/bin/gemma";
open FT,$flat;
my %ch;
while (<FT>) {
	chomp;
	next if(/^geneName/);
	my (undef,$gene,$chr,undef,$start,$end,undef)=split/\s+/,$_,7;
	$chr=~s/chr//g;
	$ch{$gene}=join("\t",$chr,$start,$end);
}
close FT;
open Pos,">$pos";
print Pos "WGT\tID\tCHR\tP0\tP1\n";
open EX,$exp;
while (<EX>) {
	chomp;
	next if(/^sample/);
	my ($id,$all)=split/\s+/,$_,2;
	my @all=split/\s+/,$all;
	open FM,"<$bfile.fam";
	open Out,">$fout/pop.fam"; 
	my $nn=0;
	while (<FM>) {
		chomp;
		my ($id1,$id2,$zr1,$zr2,$zr3,undef)=split/\s+/,$_;
		print Out "$id1\t$id2\t$zr1\t$zr2\t$zr3\t$all[$nn]\n";
		$nn++;
	}
	close FM;
	close Out;
	print Pos "$wp/$id.wgt.RDat\t$id\t$ch{$id}\n";
	`mv $fout/pop.fam pop.fam`;
	`Rscript $weight --bfile $bfile --tmp $fout/$id --out $fout/$id --models $models --PATH_gemma $gemma  --PATH_gcta $gcta`;
    `mv $fout/$id.wgt.RDat $wp`;
}
close Pos;

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

	eg: perl $Script -expression KFtpm.norm.txt -bfile pop -out ./result -weight ./weight -pos Bna.pos -refflat Brassica_napus.refflat
	
Usage:
  Options:
	"expression:s"=>\$exp,
	"bfile:s"=>\$bfile,
	"out:s"=>\$fout,
	"models:s"=>\$model,
	"weight:s"=>\$wp,
	"pos:s"=>\$pos,
	"refflat:s"=>\$flat,

USAGE
        print $usage;
        exit;
}
