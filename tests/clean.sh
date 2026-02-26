for i in $(find . -name "*.nim");
do
  bin_name=$(dirname $i)/$(basename ${i%.nim})
  if [[ -f "$bin_name" ]]; then 
    echo " !> $bin_name"
    rm "$bin_name"
  fi
done
