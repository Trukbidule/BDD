NB_COL = 10
ESC_CHAR = '#'
FILE_NAME = "./tt_test.txt"
COL_TARGET = 9
INV_RES = 1

res = []
str_res= ""

#open file and read line by line
f = open("./tt_test.txt", "r")
line = f.readline()
while line:
    
    #extract the truth table if not a comment
    if(line[0] != ESC_CHAR):
        res.append((line.split())[COL_TARGET])
    
    line = f.readline() #read next line

f.close() 

#reverse the TT is needed
if(INV_RES):
    res.reverse()

for s in res :
    str_res = str_res + str(s)
print(str_res)
