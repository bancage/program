#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use feature qw{say};

# 一个参数:文件名
# 返回文件内容，去掉所有空白符
sub get_content{
    my $file=shift;
    open my $in,"<",$file;
    $_=join "",<$in>;
    s/\s//g;
    close $in;
    return $_;
}

# 第一个参数:dna文件名
# 第二个参数:l_mrna文件名
# 返回一个数组，每个元素是a..b的一个字符串
sub get_l_mrna_range{
    my $dna_file=shift;
    my $l_mrna_file=shift;
    my $mrna_length=get_content($l_mrna_file);
    my $dna_exon_content=&get_content($dna_file);
    $dna_exon_content =~ s/<//;
    my @dna_range=map {[ split /\.\./,$_ ]} split ",",$dna_exon_content;
    my $total=0;
    my @intron_len;
    # 得到每个intron的长度
    for (my ($i,$len)=(0,scalar @dna_range);$i<$len - 1;$i++) {
    	push @intron_len,($dna_range[$i+1]->[0] - $dna_range[$i]->[1] - 1);
    }
    foreach ( @dna_range ) {
	$total += $_->[1]-$_->[0]+1;
    }
    $dna_range[0]->[0]+=$total-$mrna_length;
    my ($base,@exon_range)=(1);
    foreach (@dna_range) {
	my $len = $_->[1]-$_->[0];
	push @exon_range,"$base..@{[($base+$len)]}";
	$base+=$len+1;
    }
    return (\@exon_range,\@intron_len);
}


# 第一个参数:dna文件名
# 第二个参数:l_mrna文件名
# 返回结果字符串，调用get_l_mrna_range，然后转换成字符串返回
sub print_get_l_mrna_range{
    my $dna_file=shift;
    my $l_mrna_file=shift;
    my $output_file=shift;
    my @result=&get_l_mrna_range($dna_file,$l_mrna_file);
    my @exon_range=@{$result[0]};
    my @intron_len=@{$result[1]};
    my $i;
    open my $out,">",$output_file;
    say $out "exon number :";
    say $out scalar @exon_range;
    say $out "exon position :";
    $i=1;
    foreach(@exon_range){
	say $out "exon $i: $_";
	$i++
    }
    say $out "intron number :";
    say $out scalar @intron_len;
    say $out "intron size :";
    $i=1;
    foreach(@intron_len){
	say $out "intron $i: $_";
	$i++
    }
}
    if (@ARGV<3) {
	say "$0 : dna_file l_mrna_file output_file";
	exit(1);
    }
my $dna_file=shift;
my $l_mrna_file=shift;
my $output_file=shift;
print_get_l_mrna_range($dna_file,$l_mrna_file,$output_file);

