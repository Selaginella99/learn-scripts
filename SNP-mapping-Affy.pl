#!/usr/bin/perl -w
# 中文注释

# 在windows安装perl后，在cmd窗口 运行第四行的命令行：程序文件名，输入文件1，输入文件2，输出文件
# info_extract_fromAffy2GenAge.pl anno_chr1_22_Affy_550K.txt GenAge.txt GeneAge_anno_chr1_22_Affy_550K.txt
# 数组@AGRV是perl内置数组，用于存储从命令行中读入的文件，$ARGV[0]存储输入文件1，$ARGV[1]存储输入文件2，$ARGV[2]存储第三个文件，依次类推

open(IN1, "$ARGV[0]") || die "Unable to open $ARGV[0]: $!";   # open 函数 open (文件句柄，文件名或用于存储文件的标量)
open (IN2, "$ARGV[1]")|| die "Unable to open $ARGV[1]: $!";
open(OUTPUT, ">$ARGV[2]") || die "Unable to write to $ARGV[2]: $!"; # 输出文件跟输入文件的不同在于前面有个 “>”符号

# <IN1> 一行一行读入文件
$header = <IN1>;  # 第一次读入第一行,标头
while ($line= <IN1>){ # while循环，将每一行存储到标量$line中
    chomp $line; # chomp函数， 删除换行符“\n”;
    (undef,$snp_id,undef,undef,undef,undef,undef,undef,undef,$gene,$type)=split/\t/,$line; # split函数，以“\t"分割，将第2列，第10、11列分别存储到$snp_id,$gene 和$type
	if($gene eq ""){next;} # 如果没有基因，next：不运行以下代码，直接进入下一个循环
	$snp{$snp_id}++; # 哈希结构，%snp；将$snp_id作为关键字存储到%snp中，键值从0累加
	# @type 用于存储 3utr, coding, intron 等类型
	if(!exists($val{$type})){ # 又一个哈希表 %val，exists()函数，判断这个类型是否出现过
	    push @type,$type; #第一次出现这种类型，不存在，将这个类型存储到 @type数组中
	}
	$val{$type}++; # 从0累加这个$type的键值，下一行再出现这种类型，exists()就是yes了
	${$gene}{$type}++; # 使用“引用”，每个$gene是对应一个哈希表，每个哈希表里该基因存储出现过的类型，有多少基因，就产生了多少哈希表！！
	push @{$gene."-".$type},$snp_id; #使用“引用”，将每个snp分别存储到对应基因对应类型的数组中，所以产生了gene_number * type_number个数组！！！
}

$header = <IN2>; # 将文件二的第一行存储到$header, 注意到该标量跟第10行重名，这时$header就代替成存储文件二第一行的内容
chomp $header; # 删除换行符“\n”;
print OUTPUT "$header"; # 打印第一行
#在第一行后面紧接着打印基因的各个类型
foreach $type(sort{$a cmp $b}@type){  #sort函数，sort{$a cmp $b} 按字符串ASCII值从小到大排序，将每个类型取出存储到标量$type
    print OUTPUT "\t$type"; # 打印制表符和类型
}
print OUTPUT "\n"; #打印 换行符

while ($line= <IN2>){
    chomp $line;
    @F=split/\t/,$line; #将$line以"\t"分割，每一列依次存储到@F数组中 $F[0]存储第一列，$F[1]存储第二列，依次类推
	$gene=$F[0]; # 第一列即基因名
	print OUTPUT "$line\t"; # 打印该行
	foreach $type(sort{$a cmp $b}@type){ # 同上，对@type各类型排序，与第一行的标头对应
	    if(!exists${$gene}{$type}){} #如果这个基因的这个类型不存在，不做处理
		elsif(@{$gene."-".$type}==1){ # 如果有一个SNP 在这个基因的这个类型中
		    print OUTPUT ${$gene."-".$type}[0]; # 打印
		}
		else{    # 如果有两个以上SNP
		    print OUTPUT ${$gene."-".$type}[0]; #先打印第一个
	        foreach $snp_id (1..$#{$gene."-".$type}){ #后面的SNP
		        print OUTPUT ",",${$gene."-".$type}[$snp_id]; #先打印逗号，再打印SNP
			}	
		}
	     print OUTPUT "\t"; # 打印制表符
	}
	print OUTPUT "\n"; # 打印换行符
}
