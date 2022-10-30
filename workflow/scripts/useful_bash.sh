#!/usr/bin/env bash

function hg19_names_to_hs37d5(){
    sed s,chr,, $1 | perl -pe 's/[^\s_]+_([^\s_]+)_random/$1.1/' | tr "gl" "GL" 
}