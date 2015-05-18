for i in $(ls *.rzdat); do h2root $i; done
for i in $(ls *.hbook); do h2root $i; done
echo " "
echo " "
echo "<Marco>: Converted files (check yourself for success):"
for i in $(ls *.rzdat); do echo "$i"; done
for i in $(ls *.hbook); do echo "$i"; done
