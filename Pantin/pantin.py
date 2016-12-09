


pantin = {'00':{'nombre':'Julio Pantin Alfonso','nacio':'Caracas, Venezuela'},
          '01':{'nombre':'Federico Lewis Pantin','nacio':'Caracas, Venezuela'},
          '02':{'nombre':'Frederick Lewis Pantin Ganteaume','nacio':'Port of Spain, Trinidad'},
          '03':{'nombre':'Charles Georges John Pantin de Mouillebert','nacio':'Pointe a Pierre, Trinidad'},
          '04':{'nombre':'Lewis (IV) Pantin','nacio':'Knaresborough,  England'},
          '05':{'nombre':'Lewis (III) Pantin','nacio':'London, England'},
          '06':{'nombre':'Lewis (II) Pantin','nacio':'London, England'},
          '07':{'nombre':'Lewis (I) Armand Pantin','nacio':'London, England'},
          '08':{'nombre':'Simon Pantin','nacio': 'Rouen, France'},
          '09':{'nombre':'Isaie Pantin II','nacio':'London, England'},
          '10':{'nombre':'Isaie Pantin (I)','nacio':'Rouen, France'}}


for i in range(10):
  pantin['%02d'%i]['from']=pantin['%02d'%(i+1)]['nacio']

pantin['10']['from']= pantin['10']['nacio']


for key in map(lambda x:'%02d'%x, range(11)):
  print "------ generacion %s"%key
  print pantin[key]



















