import pandas

#rows = [1]
#for x in range(27,63):
    #rows.append(x)

#df = pandas.read_csv('standard-element-concentrations.csv', usecols = [3,4,5,6,7,8,10,12,15,18,23], skiprows=rows)

df = pandas.read_csv('standard-element-concentrations.csv')
print(df)
#print(df == 46.7400)