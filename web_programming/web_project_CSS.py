#!/usr/local/Python-3.7/bin/python
import pymysql
import os,sys
import cgi
import matplotlib
import numpy as np
import pandas as pd

import cgitb


cgitb.enable()

matplotlib.use('Agg')

os.environ['HOME'] = '/tmp'



import matplotlib.pyplot as plt
# print content-type





print("Content-type: text/html\n")

print("<html><head>")
print("<title>Enter New Order</title>")
print('''<style>
body {margin:30;padding:30}
#Products {
  font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;
  border-collapse: collapse;
  width: 50%;
}


</style>
</head>''')
print("<body>")

print("<h1>Enter New Order</h1>")

query = """SELECT id, name from Product"""
Taxonomiclist = ["species","genus","family","order","class","phylum"]
MouseType = ["Pure","Mix","All"]
DiseaseStatus = ["WT","AD","All"]


print('''<form action="https://bioed.bu.edu/cgi-bin/students_20/mkouzmin/add_order.py" method="POST" >
	Rank:<select name="Rank">
    <option value="Species">Species</option>
    <option value="Genus">Genus</option>
    <option value="Family">Family</option>
    <option value="Order">Order</option>
    <option value="Class">Class</option>
    <option value="Phylum">Phylum</option>
    </select><br />
    Mouse Type::<select name="Mouse Type">
    <option value="Pure">Pure</option>
    <option value="Mix">Mix</option>
    </select><br />
    Disease Status::<select name="Disease Status">
    <option value="AD">Alzheimer's Disease</option>
    <option value="Wild Type">Wild Type</option>
    </select><br />
	<input type="submit" value="View HeatMap">
	</form>''')

# get the form
form = cgi.FieldStorage()
if form:
    Rank = form.getvalue("Rank")
    Type = form.getvalue("Mouse Type")
    Status = form.getvalue("Disease Status")
    #print(Rank)
    #print(Type)
    #print(Status)
    query = ""

    if Type =="All" & Status == "All":
        query = """SELECT Abundance.value, mouse.name, TaxonomicRank.name FROM Abundance join mouse using mid join TaxonomicRank using tid WHERE rank ='%s';""" % (
            Rank)
    elif Type == "All":
        query = """SELECT Abundance.value, mouse.name, TaxonomicRank.name FROM Abundance join mouse using mid join TaxonomicRank using tid WHERE rank ='%s', WTvsAD = '%s';""" % (
        Rank, Status)
    elif Status == "All":
        query = """SELECT Abundance.value, mouse.name, TaxonomicRank.name FROM Abundance join mouse using mid join TaxonomicRank using tid WHERE rank ='%s', PurevsMix = '%s';""" % (
            Rank, Type)
    else:
        query = """SELECT Abundance.value, mouse.name, TaxonomicRank.name FROM Abundance join mouse using mid join TaxonomicRank using tid WHERE rank ='%s', PurevsMix = '%s', WTvsAD = '%s';"""% (Rank,Type,Status)


    # print(query)
    Abundance = []
    MName = []
    TName = []

    try:
        connection = pymysql.connect(host="bioed.bu.edu", user="mkouzmin", password="mkouzmin", db="mkouzmin", port=4253)
        df = pd.read_sql(query,connection)
        table = df.pivot(index='TaxonomicRank.name', columns='mouse.name', values='Abundance.value')
        fig, ax = plt.subplots()
        ax.pcolor(table.values, cmap=plt.get_cmap('jet'),
                  vmin=df['min_value'].min(), vmax=df['min_value'].max())
        ax.set_xticks(np.arange(table.shape[1] + 1) + 0.5, minor=False)
        ax.set_xticklabels(table.columns, minor=False)
        ax.set_yticks(np.arange(table.shape[0] + 1) + 0.5, minor=False)
        ax.set_yticklabels(table.index, minor=False)
        ax.set_xlim(0, table.shape[1])
        ax.set_ylim(0, table.shape[0])

        print("Content-type: image/png")
        print()
        plt.savefig(sys.stdout, format='png')
    except Exception as mysqlError:
        print("<p><font color=red><b>Error</b> while executing query</font></p>")


print("</body></html>")
