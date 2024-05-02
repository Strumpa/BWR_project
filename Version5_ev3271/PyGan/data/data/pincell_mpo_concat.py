#
# Concatenation of MPO files
#
import lifo
import cle2000
#
# Call a Cle2000 procedure to perform concatenation
ipLifo=lifo.new()
ipLifo.push("MyMpo1.hdf")
ipLifo.push("MyMpo2.hdf")
ipLifo.push("MyMpoCat.hdf")
ConcatProc = cle2000.new('MPOconcat',ipLifo,1)
ConcatProc.exec()
print("test pincell_mpo_concat completed")
