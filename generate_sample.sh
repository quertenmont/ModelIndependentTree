executable=`echo $0 | sed 's/.sh/.C/'` #assume the .sh and .C file have the same name
args=''
for var in "$@"
do
    if [ "$args" ]; then args=$args, ; fi
    args=$args\"$var\"
done

echo $args

root -l -b << EOF 
  TString makeshared(gSystem->GetMakeSharedLib());
  makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated ");
  makeshared.ReplaceAll("-Wshadow ", " -std=c++0x -D__USE_XOPEN2K8 ");
  gSystem->SetMakeSharedLib(makeshared);
  .x $executable+($args)
EOF

