#!/bin/bash


infile=$1
toward=$2

echo $infile

if  [ $toward = "num" ];then
sed "s/NC_053212.1/1/" $infile |
    sed "s/NC_053213.1/2/" |
    sed "s/NC_053214.1/3/" |
    sed "s/NC_053215.1/4/"|
    sed "s/NC_053216.1/5/"|
    sed "s/NC_053217.1/6/"|
    sed "s/NC_053218.1/7/"|
    sed "s/NC_053219.1/8/"|
    sed "s/NC_053220.1/9/"|
    sed "s/NC_053221.1/10/"|
    sed "s/NC_053222.1/11/" |
    sed "s/NC_053223.1/12/"|
    sed "s/NC_053224.1/13/"|
    sed "s/NC_053225.1/14/"|
    sed "s/NC_053226.1/15/"|
    sed "s/NC_053227.1/16/" |
    sed "s/NC_053228.1/17/" |
    sed "s/NC_053229.1/18/" |
    sed "s/NC_053230.1/19/" |
    sed "s/NC_053231.1/20/"|
    sed "s/NC_053232.1/21/"|
    sed "s/NC_053233.1/Y/"|
    sed "s/NC_041244.1/M/"|
    sed "s/NW_[0-9\.]*/Un/"  >${infile}_num

 else
 sed "s/NC_053212.1/chrI/" $infile |
    sed "s/NC_053213.1/chrII/" |
    sed "s/NC_053214.1/chrIII/" |
    sed "s/NC_053215.1/chrIV/"|
    sed "s/NC_053216.1/chrV/"|
    sed "s/NC_053217.1/chrVI/"|
    sed "s/NC_053218.1/chrVII/"|
    sed "s/NC_053219.1/chrVIII/"|
    sed "s/NC_053220.1/chrIX/"|
    sed "s/NC_053221.1/chrX/"|
    sed "s/NC_053222.1/chrXI/" |
    sed "s/NC_053223.1/chrXII/"|
    sed "s/NC_053224.1/chrXIII/"|
    sed "s/NC_053225.1/chrXIV/"|
    sed "s/NC_053226.1/chrXV/"|
    sed "s/NC_053227.1/chrXVI/" |
    sed "s/NC_053228.1/chrXVII/" |
    sed "s/NC_053229.1/chrXVIII/" |
    sed "s/NC_053230.1/chrXIX/" |
    sed "s/NC_053231.1/chrXX/"|
    sed "s/NC_053232.1/chrXXI/"|
    sed "s/NC_053233.1/chrY/"|
    sed "s/NC_041244.1/chrM/">${infile}_chr
    #sed "s/NW_[0-9\.]*/chrUn/"  >${infile}_chr
fi
