ccol=~/src/prog/cCOLUMBUS 
echo molint
${ccol}/molint < int.in > int.out
echo cdip
${ccol}/cdip < int.in > cdip.out
echo pkfl
${ccol}/pkfl < dum.in > pkfl.out
echo scf
${ccol}/cscf < scf.in > scf.out
echo "CI(Groud)"
${ccol}/molci < ci1.in > ci1.out
mv CSF CSFG
mv CIVEC CIVECG
echo "CI(EE)"
${ccol}/molci < ci2.in > ci2.out
