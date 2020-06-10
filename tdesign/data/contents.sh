#!/bin/bash
for entry in *
do
  echo "\"../data/$entry\","
done


echo ""
for entry in *
do
  echo "${entry##*.},"
done
