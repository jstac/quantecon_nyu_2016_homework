import matplotlib.pyplot as plt
import requests

data_file=open('company_list.csv').readlines()
bad_chars = '+$"%\n' #Eliminate these characters
header= data_file.pop(0) # Remove the header
mkc=[]  #initialization
symb=[]
chg=[]
counter=0 # Just to see where the downloading is at
for line in data_file:
    symbol, name, marketcap = line.split(',') # Split the lines
    for c in bad_chars: # Clean the lines
       marketcap=marketcap.replace(c,"")
       symbol=symbol.replace(c,"")
    marketcap=marketcap.replace("n/a","0") # convention
    marketcap=marketcap.replace("Ltd.","0") # convention
    marketcap=marketcap.replace("Inc.","0") # convention
    if(marketcap.find("M")==-1): #If millions, multiply value by million
        if(marketcap.find("B")==-1):
            marketcap=float(marketcap)
            if(marketcap<1000): 
                marketcap=marketcap*1000000 #Assume that it is always millions
        else:
            marketcap=marketcap.replace("B","")
            marketcap=float(marketcap)
            marketcap=1000000000*marketcap    
    else:
        marketcap=marketcap.replace("M","")
        marketcap=float(marketcap)
        marketcap=1000000*marketcap
    marketcap=marketcap/1000000   # Everything is in millions
    if(marketcap!=0): # Ignore data otherwise
        mkc.append(marketcap)
        symb.append(symbol)
        counter+=1
        link='https://finance.yahoo.com/d/quotes.csv?s=%s&f=np2' %symbol
        data=requests.get(link)
        print(counter)
        stringchange=data.text
        print( stringchange)
        change=stringchange.rpartition(',')[-1] # Tricky -save last part
        for c in bad_chars:
            change=change.replace(c,"")
        change=float(change)
        chg.append(change)
plt.scatter(mkc,chg) #Shows that smaller firms are more volatile
plt.show()
