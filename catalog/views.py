from django.shortcuts import render
from .models import Experiment, Geneexpression, Geneset, Patient
import json
import pandas as pd
from django.db import connection
import numpy as np
from dal import autocomplete
import plotly.graph_objs as go
import plotly.offline as opy
from plotly.subplots import make_subplots

#Initializing website
geneNames = Geneexpression.objects.all()
newGeneNames = set()
for k in geneNames:
    newGeneNames.add(k.genename)
newGeneNames = list(newGeneNames)
newGeneNames.sort()

allSamples = ['1672MH8C','1672MX24C','1672MY25C', '1672MI9T','1672MK11T','1672MM13T', '2011AC3C',
              '2011AF6C','2011AG7C','2011AH8C','2011AA1T','2011AB2T', '2011AD4T', '2011AE5T',
              '272AG7C','272AP16C','272AW23C','272AZ26C','272AB2T','272ACC29T', '272ADD30T',
              '272AS19T', '2012ADD30C','2012AF6C','2012AR18C','2012AX24C','2012AD4T','2012AJ10T',
              '2012AW23T','2012AY25T']

rmSampleList = []
currentSamples = ['1672MH8C','1672MX24C','1672MY25C', '1672MI9T','1672MK11T','1672MM13T', '2011AC3C',
              '2011AF6C','2011AG7C','2011AH8C','2011AA1T','2011AB2T', '2011AD4T', '2011AE5T',
              '272AG7C','272AP16C','272AW23C','272AZ26C','272AB2T','272ACC29T', '272ADD30T',
              '272AS19T', '2012ADD30C','2012AF6C','2012AR18C','2012AX24C','2012AD4T','2012AJ10T',
              '2012AW23T','2012AY25T']

def index(request):
    """
    View function for home page of site.
    """
    return render(
        request,
        'index.html',
        {'searchGeneNames': json.dumps(newGeneNames)}
    )

#Handling plotting request.
def plot(request):
    if request.method == 'POST':
        if 'plot' in request.POST:
            query = request.POST.get('plotQ')
            global plotType
            plotType = request.POST.get('plotType')
            file = request.FILES.get('fileData')
            global ar
            ar = request.POST.get('groupBy')
            global geneList
            geneList = list()
            notMatchGenes = list()
            global currentSamples
            currentSamples = ['1672MH8C','1672MX24C','1672MY25C', '1672MI9T','1672MK11T','1672MM13T', '2011AC3C',
              '2011AF6C','2011AG7C','2011AH8C','2011AA1T','2011AB2T', '2011AD4T', '2011AE5T',
              '272AG7C','272AP16C','272AW23C','272AZ26C','272AB2T','272ACC29T', '272ADD30T',
              '272AS19T', '2012ADD30C','2012AF6C','2012AR18C','2012AX24C','2012AD4T','2012AJ10T',
              '2012AW23T','2012AY25T']
            #If users type gene names by themselves
            if query:
                query = query.splitlines()
                text = list()
                for i in query:
                    i = i.strip()
                    text.append(i)
            #If users upload a file
            elif file:
                if file.name.endswith('.txt') or file.name.endswith('.csv'):
                    text = pd.read_csv(file, header=None)
                    text = text.iloc[:, 0]
                elif file.name.endswith('.xlsx'):
                    text = pd.read_excel(file, header=None)
                    text = text.iloc[:, 0]
                else:
                    error = 'File is not .txt, .csv or .xlsx type'
                    return render(request, 'index.html', {'errorFile': error})
            #Processing typed gene or uploaded file. Both of them are called text.
            for i in text:
                gene = Geneexpression.objects.filter(genename__iexact=i)
                if len(gene) != 0:
                    # Check whether the gene is in the database
                    for j in gene:
                        geneList.append(j.genename)
                else:
                    notMatchGenes.append(i)

            #convertListToString is a function to convert a list to a string
            notMatchGenes = convertListToString(notMatchGenes)
            info = 'Not match gene(s): ' + notMatchGenes
            geneList = list(set(geneList))
            if len(geneList) != 0:
                if plotType == 'boxplotT':
                    src = boxplot(geneList, ar, currentSamples)
                elif plotType == 'heatmapT':
                    if ar == 'noAr':
                        src = heatmap(geneList, currentSamples)
                    else:
                        src = heatmapAr(geneList, currentSamples, ar)
                if len(notMatchGenes) != 0:
                    context = {'plot': src, 'notMatch': info, 'geneList': json.dumps(geneList),
                               'searchGeneNames': json.dumps(newGeneNames),
                               'currentSamples':json.dumps(allSamples), 'plotType':plotType}
                else:
                    context = {'plot': src, 'geneList': json.dumps(geneList),
                               'searchGeneNames': json.dumps(newGeneNames),
                               'currentSamples':json.dumps(allSamples), 'plotType':plotType}
                return render(request, 'plot.html', context)
            else:
                return render(request, 'plot.html', {'notMatch': info})

#Request is plotting boxplot. Initial processing data.
def boxplot(gList, ar, sp):
    if ar == 'noAr':
        sql = (' SELECT g.geneName, g.geneExpressionValues, e.exIdName '
           ' FROM geneExpression g, experiment e '
           ' where e.exId = g.exId and g.geneName in (%s)'
           % (','.join(['%s'] * len(gList))))
    else:
        sql = ('SELECT p.arType, g.geneName, g.geneExpressionValues, e.type, e.exIdName '
               'FROM patient p, experiment e, geneExpression g '
               'where p.patientId = e.patientId and e.exId = g.exId '
               'and g.geneName in (%s)'
               % (','.join(['%s'] * len(gList))))

    # Create dataframe
    cur = connection.cursor()
    cur.execute(sql, gList)

    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    global boxdf
    boxdf = pd.DataFrame(rows, columns=names)

    # Convert gene expression values to logarithm
    boxdf["geneExpressionValues"] += 0.1
    boxdf["geneExpressionValues"] = np.log2(boxdf["geneExpressionValues"])
    if len(sp) != 30:
        boxdf = boxdf.loc[boxdf['exIdName'].isin(sp)]
    div = plotBoxplot(boxdf, ar)
    return div

#Plotting boxplot
def plotBoxplot(df, ar):
    x1 = df['geneName'].tolist()
    # Ploting
    if ar == 'noAr':
        y1 = df['geneExpressionValues'].tolist()
        trace1 = go.Box(x=x1, y=y1)
        data = go.Data([trace1])
        figure = go.Figure(data=data)
        figure.update_layout(xaxis=dict(title_text="Gene Name", title_font={"size": 17}, title_standoff=25),
                             yaxis=dict(title_text="Gene Expression Values (log2)", title_standoff=25),
                             hoverlabel=dict(align='left'))
    elif ar == 'withAr':
        newArPdf = df[boxdf['arType'] == 'AR+']
        newArNdf = df[boxdf['arType'] == 'AR-']
        yArP = newArPdf['geneExpressionValues'].tolist()
        yArN = newArNdf['geneExpressionValues'].tolist()
        trace1 = go.Box(x=x1, y=yArP, name='AR+')
        trace2 = go.Box(x=x1, y=yArN, name='AR-')
        data = go.Data([trace1, trace2])

    elif ar == 'withSample':
        newSampleCdf = df[boxdf['type'] == 'control']
        newSampleTdf = df[boxdf['type'] == 'treatment']
        ySampleC = newSampleCdf['geneExpressionValues'].tolist()
        ySampleT = newSampleTdf['geneExpressionValues'].tolist()
        trace1 = go.Box(x=x1, y=ySampleC, name='Control')
        trace2 = go.Box(x=x1, y=ySampleT, name='Treatment')
        data = go.Data([trace1, trace2])

    layout = go.Layout(boxmode='group')
    figure = go.Figure(data=data, layout=layout)
    figure.update_layout(xaxis=dict(title_text="Gene Name", title_font={"size": 17}, title_standoff=25),
                         yaxis=dict(title_text="Gene Expression Values (log2)", title_standoff=25))

    div = opy.plot(figure, auto_open=False, output_type='div')
    return div

#Request is plotting heatmap. Initial selecting and processing data.
def heatmap(gList, sp):
    sqlNoAr = "SELECT g.geneName, g.geneExpressionValues, e.exIdName FROM patient p, experiment e, geneExpression g where p.patientId = e.patientId and e.exId = g.exId and g.geneName in (%s)" % (
        ','.join(['%s'] * len(gList)))

    # Create a dataframe without AR type
    cur = connection.cursor()
    cur.execute(sqlNoAr, gList)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    # Create dataframe without AR types for multiple genes

    global noAr
    noAr = pd.DataFrame(rows, columns=names).set_index('exIdName')
    noAr["geneExpressionValues"] += 0.1
    noAr["geneExpressionValues"] = np.log2(noAr["geneExpressionValues"])
    if len(sp) != 30:
        noAr = noAr[noAr.index.isin(sp)]
    div = plotHeatmapNoAr(gList, noAr)
    return div

#Plotting heatmap without AR and sample types.
def plotHeatmapNoAr(geneList,samplesDf):
    fig = go.Figure()
    #If there are more than one genes
    if len(geneList) != 1:
        fig.add_trace(go.Heatmap(
            z=samplesDf["geneExpressionValues"],
            x=samplesDf["geneName"],
            y=samplesDf.index,
            name='',
            #coloraxis="coloraxis",
            colorscale=[[0.0, "rgb(49,54,149)"],
                        [0.1111111111111111, "rgb(69,117,180)"],
                        [0.2222222222222222, "rgb(116,173,209)"],
                        [0.3333333333333333, "rgb(171,217,233)"],
                        [0.4444444444444444, "rgb(224,243,248)"],
                        [0.5555555555555556, "rgb(254,224,144)"],
                        [0.6666666666666666, "rgb(253,174,97)"],
                        [0.7777777777777778, "rgb(244,109,67)"],
                        [0.8888888888888888, "rgb(215,48,39)"],
                        [1.0, "rgb(165,0,38)"]],
            colorbar=dict(title='Log2 Scale'),
            hovertemplate='Gene Name: %{x} <br>Sample Name: %{y}<br>Gene Expression Value: %{z}'))
        fig.update_layout(autosize=False, width=700, height=800)
        fig.data[0].update(zmin=-4, zmax=21)
    #If there is only one gene
    else:
        fig.add_trace(go.Heatmap(
            z=samplesDf["geneExpressionValues"],
            x=samplesDf.index,
            y=samplesDf["geneName"],
            name='',
            colorscale=[[0.0, "rgb(49,54,149)"],
                        [0.1111111111111111, "rgb(69,117,180)"],
                        [0.2222222222222222, "rgb(116,173,209)"],
                        [0.3333333333333333, "rgb(171,217,233)"],
                        [0.4444444444444444, "rgb(224,243,248)"],
                        [0.5555555555555556, "rgb(254,224,144)"],
                        [0.6666666666666666, "rgb(253,174,97)"],
                        [0.7777777777777778, "rgb(244,109,67)"],
                        [0.8888888888888888, "rgb(215,48,39)"],
                        [1.0, "rgb(165,0,38)"]],
            colorbar=dict(title='Log2 Scale'),
            hovertemplate='Gene Name: %{y} <br>Sample Name: %{x}<br>Gene Expression Value: %{z}'))
        fig.update_layout(autosize=False, width=800, height=280)
        fig.data[0].update(zmin=-4, zmax=21)
    div = opy.plot(fig, auto_open=False, output_type='div')
    return div

#Initial processing data for heatmap with AR or sample types, or both of them.
def heatmapAr(gList, sp, ar):
    sqlNoAr = "SELECT g.geneName, g.geneExpressionValues, e.exIdName FROM patient p, experiment e, geneExpression g where p.patientId = e.patientId and e.exId = g.exId and g.geneName in (%s)" % (','.join(['%s'] * len(gList)))
    sqlAr = 'SELECT e.exIdName, e.type, p.arType FROM patient p, experiment e where p.patientId = e.patientId'

    # Create a dataframe without AR type
    cur = connection.cursor()
    cur.execute(sqlNoAr, gList)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    # Create dataframe without AR types for multiple genes
    global dfSamples
    dfSamples = pd.DataFrame(rows, columns=names).set_index('exIdName')
    dfSamples["geneExpressionValues"] = np.log2(dfSamples["geneExpressionValues"] + 0.1)

    global dfAr
    dfAr = pd.read_sql(sqlAr, connection, index_col="exIdName")
    #If sample size is not 30, filter removed or added sample from the dataframe
    if len(sp) != 30:
        dfSamples = dfSamples[dfSamples.index.isin(sp)]
        dfAr = dfAr[dfAr.index.isin(sp)]
    div = plotHeatmapAr(gList, dfSamples, dfAr, ar)
    return div

#Modifying data and plotting heatmap with AR or sample types or both of them.
def plotHeatmapAr(gList,samplesDf, arDf, ar):
    if len(gList) != 1:
        if 'withAr' or 'allTypes' in ar:
            newArList = []
            newArStrList = []
            for l in arDf["arType"].to_list():
                sub = l.split(',')
                newArList.append(sub)

            for n in newArList:
                while 'AR+' in n:
                    n[n.index('AR+')] = 20
                while 'AR-' in n:
                    n[n.index('AR-')] = -3

            for l in arDf["arType"].to_list():
                sub = l.split(',')
                newArStrList.append(sub)
        if 'withSample' or 'allTypes' in ar:
            typeList = []
            typeStrList = []
            for l in arDf["type"].to_list():
                t = l.split(',')
                typeList.append(t)

            for n in typeList:
                while 'control' in n:
                    n[n.index('control')] = 1
                while 'treatment' in n:
                    n[n.index('treatment')] = 14

            for l in arDf["type"].to_list():
                b = l.split(',')
                typeStrList.append(b)

        if ar == 'allTypes':
            fig = make_subplots(1, 3, column_widths=[0.08, 0.08, 0.84], shared_yaxes=True)
        else:
            fig = make_subplots(1, 2, column_widths=[0.1, 0.9], shared_yaxes=True)

        if ar == 'withAr':
            fig.add_trace(go.Heatmap(
                                z=newArList,
                                x=['AR Type'],
                                y=arDf.index,
                                text=newArStrList,
                                hovertemplate='Sample Name: %{y}<br>AR Type: %{text}',
                                colorscale='Cividis',
                                name=''), 1, 1)
            fig.add_trace(go.Heatmap(
                                z=samplesDf["geneExpressionValues"],
                                x=samplesDf["geneName"],
                                y=samplesDf.index,
                                hovertemplate='Gene Name: %{x} <br>Sample Name: %{y}<br>Gene Expression Value: %{z}',
                                #Setting the heatmap color
                                colorscale=[[0.0, "rgb(49,54,149)"],
                                            [0.1111111111111111, "rgb(69,117,180)"],
                                            [0.2222222222222222, "rgb(116,173,209)"],
                                            [0.3333333333333333, "rgb(171,217,233)"],
                                            [0.4444444444444444, "rgb(224,243,248)"],
                                            [0.5555555555555556, "rgb(254,224,144)"],
                                            [0.6666666666666666, "rgb(253,174,97)"],
                                            [0.7777777777777778, "rgb(244,109,67)"],
                                            [0.8888888888888888, "rgb(215,48,39)"],
                                            [1.0, "rgb(165,0,38)"]],
                                colorbar=dict(title='Log2 Scale'),
                                name=''), 1, 2)
            fig.update_layout(width=700, height=700, autosize=False)
            #Hiding the colorscale of AR types
            fig.data[0].showscale = False
            #Fixing the colorscale range
            fig.data[0].update(zmin=-4, zmax=21)
            fig.data[1].update(zmin=-4, zmax=21)
        elif ar == 'withSample':
            fig.add_trace(go.Heatmap(
                                z=typeList,
                                x=['Sample Type'],
                                y=arDf.index,
                                text=typeStrList,
                                hovertemplate='Sample Name: %{y}<br>Sample Type: %{text}',
                                colorscale="Blackbody",
                                name=''), 1, 1)
            fig.add_trace(go.Heatmap(
                                z=samplesDf["geneExpressionValues"],
                                x=samplesDf["geneName"],
                                y=samplesDf.index,
                                hovertemplate='Gene Name: %{x} <br>Sample Name: %{y}<br>Gene Expression Value: %{z}',
                                colorscale=[[0.0, "rgb(49,54,149)"],
                                            [0.1111111111111111, "rgb(69,117,180)"],
                                            [0.2222222222222222, "rgb(116,173,209)"],
                                            [0.3333333333333333, "rgb(171,217,233)"],
                                            [0.4444444444444444, "rgb(224,243,248)"],
                                            [0.5555555555555556, "rgb(254,224,144)"],
                                            [0.6666666666666666, "rgb(253,174,97)"],
                                            [0.7777777777777778, "rgb(244,109,67)"],
                                            [0.8888888888888888, "rgb(215,48,39)"],
                                            [1.0, "rgb(165,0,38)"]],
                                colorbar=dict(title='Log2 Scale'),
                                name=''), 1, 2)
            fig.update_layout(width=700, height=700, autosize=False)
            fig.data[0].showscale = False
            fig.data[0].update(zmin=1, zmax=14)
            fig.data[1].update(zmin=-4, zmax=21)
        elif ar == 'allTypes':
            fig.add_trace(go.Heatmap(
                                    z=newArList,
                                    x=['AR Type'],
                                    y=arDf.index,
                                    text=newArStrList,
                                    hovertemplate='Sample Name: %{y}<br>AR Type: %{text}',
                                    colorscale="Cividis",
                                    name=''), 1, 1)
            fig.add_trace(go.Heatmap(
                                    z=typeList,
                                    x=['Sample Type'],
                                    y=arDf.index,
                                    text=typeStrList,
                                    hovertemplate='Sample Name: %{y}<br>Sample Type: %{text}',
                                    colorscale='Blackbody',
                                    name=''), 1, 2)
            fig.add_trace(go.Heatmap(
                                z=samplesDf["geneExpressionValues"],
                                x=samplesDf["geneName"],
                                y=samplesDf.index,
                                hovertemplate='Gene Name: %{x} <br>Sample Name: %{y}<br>Gene Expression Value: %{z}',
                                colorscale=[[0.0, "rgb(49,54,149)"],
                                            [0.1111111111111111, "rgb(69,117,180)"],
                                            [0.2222222222222222, "rgb(116,173,209)"],
                                            [0.3333333333333333, "rgb(171,217,233)"],
                                            [0.4444444444444444, "rgb(224,243,248)"],
                                            [0.5555555555555556, "rgb(254,224,144)"],
                                            [0.6666666666666666, "rgb(253,174,97)"],
                                            [0.7777777777777778, "rgb(244,109,67)"],
                                            [0.8888888888888888, "rgb(215,48,39)"],
                                            [1.0, "rgb(165,0,38)"]],
                                colorbar=dict(title='Log2 Scale'), name=''), 1, 3)
            fig.update_layout(width=800, height=700, autosize=False)
            #Not showing colorbars of AR and sample types
            fig.data[0].showscale = False
            fig.data[1].showscale = False
            #Setting gene expression values range so the color scale could be fixed
            fig.data[0].update(zmin=-4, zmax=21)
            fig.data[1].update(zmin=1, zmax=14)
            fig.data[2].update(zmin=-4, zmax=21)
    #If there is only one gene in the plotting list.
    else:
        if 'withAr' or 'allTypes' in ar:
            newArList = arDf["arType"].to_list()
            modifiedArList = []
            for i in newArList:
                if i == 'AR+':
                    newArList[newArList.index('AR+')] = 20
                elif i == 'AR-':
                    newArList[newArList.index('AR-')] = -3
            modifiedArList.append(newArList)
        if 'withSample' or 'allTypes' in ar:
            arTextList = arDf["arType"].to_list()
            modifiedArText = []
            modifiedArText.append(arTextList)

            typeList = arDf["type"].to_list()
            modifiedSampleTypeList = []
            for i in typeList:
                if i == 'control':
                    typeList[typeList.index('control')] = 1
                elif i == 'treatment':
                    typeList[typeList.index('treatment')] = 14
            modifiedSampleTypeList.append(typeList)

            sampleTypeTextList = arDf["type"].to_list()
            modifiedSampleText = []
            modifiedSampleText.append(sampleTypeTextList)

        #Setting subplot framework.
        if ar == 'allTypes':
            fig = make_subplots(3, 1, row_heights=[0.28, 0.28, 0.44], shared_xaxes=True)
        else:
            fig = make_subplots(2, 1, row_heights=[0.4, 0.6], shared_xaxes=True)

        if ar == 'withAr':
            fig.add_trace(go.Heatmap(
                z=modifiedArList,
                x=samplesDf.index,
                y=['AR Type'],
                text=modifiedArText,
                hovertemplate='Sample Name: %{x}<br>AR Type: %{text}',
                colorscale="Cividis",
                name=''), 1, 1)
            fig.add_trace(go.Heatmap(z=samplesDf["geneExpressionValues"],
                                     x=samplesDf.index,
                                     y=samplesDf["geneName"],
                                     colorscale=[[0.0, "rgb(49,54,149)"],
                                                 [0.1111111111111111, "rgb(69,117,180)"],
                                                 [0.2222222222222222, "rgb(116,173,209)"],
                                                 [0.3333333333333333, "rgb(171,217,233)"],
                                                 [0.4444444444444444, "rgb(224,243,248)"],
                                                 [0.5555555555555556, "rgb(254,224,144)"],
                                                 [0.6666666666666666, "rgb(253,174,97)"],
                                                 [0.7777777777777778, "rgb(244,109,67)"],
                                                 [0.8888888888888888, "rgb(215,48,39)"],
                                                 [1.0, "rgb(165,0,38)"]],
                                     colorbar=dict(title='Log2 Scale'),
                                     hovertemplate='Gene Name: %{y} <br>Sample Name: %{x}<br>Gene Expression Value: %{z}',
                                     name=''), 2, 1)
            fig.update_layout(width=700, height=280, autosize=False)
            #Hiding the colorscale of AR type
            fig.data[0].showscale = False
            #Adjusting and fixing the color
            fig.data[0].update(zmin=-4, zmax=21)
            fig.data[1].update(zmin=-4, zmax=21)
        elif ar == 'withSample':
            fig.add_trace(go.Heatmap(z=modifiedSampleTypeList,
                                    x=samplesDf.index,
                                    y=['Sample Type'],
                                    text=modifiedSampleText,
                                    hovertemplate='Sample Name: %{x}<br>AR Type: %{text}',
                                    colorscale='Blackbody',
                                    name=''), 1, 1)
            fig.add_trace(go.Heatmap(z=samplesDf["geneExpressionValues"],
                                    x=samplesDf.index,
                                    y=samplesDf["geneName"],
                                    colorscale=[[0.0, "rgb(49,54,149)"],
                                                [0.1111111111111111, "rgb(69,117,180)"],
                                                [0.2222222222222222, "rgb(116,173,209)"],
                                                [0.3333333333333333, "rgb(171,217,233)"],
                                                [0.4444444444444444, "rgb(224,243,248)"],
                                                [0.5555555555555556, "rgb(254,224,144)"],
                                                [0.6666666666666666, "rgb(253,174,97)"],
                                                [0.7777777777777778, "rgb(244,109,67)"],
                                                [0.8888888888888888, "rgb(215,48,39)"],
                                                [1.0, "rgb(165,0,38)"]],
                                    colorbar=dict(title='Log2 Scale'),
                                    hovertemplate='Gene Name: %{y} <br>Sample Name: %{x}<br>Gene Expression Value: %{z}',
                                    name=''), 2, 1)
            fig.update_layout(width=700, height=280, autosize=False)
            #Hiding the colorscale of sample types
            fig.data[0].showscale = False
            # Adjusting and fixing color range.
            fig.data[0].update(zmin=1, zmax=14)
            fig.data[1].update(zmin=-4, zmax=21)
        else:
            fig.add_trace(go.Heatmap(z=modifiedArList,
                                    x=samplesDf.index,
                                    y=['AR Type'],
                                    text=modifiedArText,
                                    hovertemplate='Sample Name: %{x}<br>AR Type: %{text}',
                                    colorscale="Cividis",
                                    name=''), 1, 1)
            fig.add_trace(go.Heatmap(z=modifiedSampleTypeList,
                                    x=samplesDf.index,
                                    y=['Sample Type'],
                                    text=modifiedSampleText,
                                    hovertemplate='Sample Name: %{x}<br>AR Type: %{text}',
                                    colorscale="Blackbody",
                                    name=''), 2, 1)
            fig.add_trace(go.Heatmap(z=samplesDf["geneExpressionValues"],
                                    x=samplesDf.index,
                                    y=samplesDf["geneName"],
                                    colorscale=[[0.0, "rgb(49,54,149)"],
                                                [0.1111111111111111, "rgb(69,117,180)"],
                                                [0.2222222222222222, "rgb(116,173,209)"],
                                                [0.3333333333333333, "rgb(171,217,233)"],
                                                [0.4444444444444444, "rgb(224,243,248)"],
                                                [0.5555555555555556, "rgb(254,224,144)"],
                                                [0.6666666666666666, "rgb(253,174,97)"],
                                                [0.7777777777777778, "rgb(244,109,67)"],
                                                [0.8888888888888888, "rgb(215,48,39)"],
                                                [1.0, "rgb(165,0,38)"]],
                                    colorbar=dict(title='Log2 Scale'),
                                    hovertemplate='Gene Name: %{y} <br>Sample Name: %{x}<br>Gene Expression Value: %{z}',
                                    name=''), 3, 1)

            fig.update_layout(width=700, height=350,autosize=False)
            fig.data[0].showscale = False
            fig.data[1].showscale = False
            #Adjusting and fixing color range.
            fig.data[0].update(zmin=-4, zmax=21)
            fig.data[1].update(zmin=1, zmax=14)
            fig.data[2].update(zmin=-4, zmax=21)
    div = opy.plot(fig, auto_open=False, output_type='div')
    return div

#Add a gene to the plot
def addGeneToPlot(request):
    if request.method == 'POST':
        addGene = request.POST.get('addGene')
        if addGene not in geneList:
            geneList.append(addGene)
        # Checking the samples that are not in the plots
        rmSampleList = list(set(currentSamples) ^ set(allSamples))
        if plotType == 'boxplotT':
            src = boxplot(geneList, ar, currentSamples)
        elif plotType == 'heatmapT':
            if ar == 'noAr':
                src = heatmap(geneList, currentSamples)
            else:
                src = heatmapAr(geneList, currentSamples, ar)
        context = {'plot': src, 'geneList': json.dumps(geneList), 'searchGeneNames': json.dumps(newGeneNames),
                   'currentSamples': json.dumps(currentSamples), 'addableSamples':json.dumps(rmSampleList), 'plotType':plotType}
        return render(request, 'plot.html', context)

#Remove a gene from the plot
def removeGeneFromPlot(request):
    if request.method == 'POST':
        if 'removeGeneButton' in request.POST:
            rmGene = request.POST.get('removeGene')
            #Checking whether the removed gene is in the plot. If so, remove it from the current plot.
            if rmGene in geneList:
                geneList.remove(rmGene)
            #Checking the samples that are not in the plots
            rmSampleList = list(set(currentSamples) ^ set(allSamples))
            if len(geneList) != 0:
                if plotType == 'boxplotT':
                    src = boxplot(geneList, ar, currentSamples)
                elif plotType == 'heatmapT':
                    if ar == 'noAr':
                        src = heatmap(geneList, currentSamples)
                    else:
                        src = heatmapAr(geneList, currentSamples, ar)
                context = {'plot': src, 'geneList': json.dumps(geneList), 'searchGeneNames': json.dumps(newGeneNames),
                           'currentSamples':json.dumps(currentSamples), 'addableSamples':json.dumps(rmSampleList), 'plotType':plotType}
                return render(request, 'plot.html', context)
            else:
                info = "No gene in the list"
                return render(request, 'plot.html', {'noGene': info, 'searchGeneNames': json.dumps(newGeneNames),
                           'currentSamples':json.dumps(currentSamples)})

#Remove a sample from the plot.
def removeSampleFromPlot(request):
    if request.method == 'POST':
        rmSample = request.POST.get('removeSample')
        #Saving the removed sample to a list so the removed samples could
        #appear at the add sample searching bar.
        if rmSample in currentSamples:
            currentSamples.remove(rmSample)
        rmSampleList = list(set(currentSamples) ^ set(allSamples))
        if len(geneList) != 0 and len(currentSamples) != 0:
            if plotType == 'boxplotT':
                src = boxplot(geneList, ar, currentSamples)
            elif plotType == 'heatmapT':
                if ar == 'noAr':
                    src = heatmap(geneList, currentSamples)
                else:
                    src = heatmapAr(geneList, currentSamples, ar)
            context = {'plot': src, 'geneList': json.dumps(geneList), 'searchGeneNames': json.dumps(newGeneNames),
                       'currentSamples':json.dumps(currentSamples), 'addableSamples':json.dumps(rmSampleList), 'plotType':plotType}
            return render(request, 'plot.html', context)
        else:
            if len(geneList) != 0 and len(currentSamples) ==0:
                info = "No sample in the list"
            elif len(geneList) == 0 and len(currentSamples) !=0:
                info = "No gene in the list"
            else:
                info = "No gene and sample in the list"
            return render(request, 'plot.html', {'noGene': info, 'geneList': json.dumps(geneList),
                                                 'searchGeneNames': json.dumps(newGeneNames),
                                                 'currentSamples':json.dumps(currentSamples),
                                                 'addableSamples':json.dumps(allSamples)})

#Add a sample to the plot.
def addSampleToPlot(request):
    if request.method == 'POST':
        addSample = request.POST.get('addSample')
        if addSample in allSamples and addSample not in currentSamples:
            currentSamples.append(addSample)
        rmSampleList = list(set(currentSamples) ^ set(allSamples))
        if len(geneList) != 0 and len(currentSamples) != 0:
            if plotType == 'boxplotT':
                src = boxplot(geneList, ar, currentSamples)
            elif plotType == 'heatmapT':
                if ar == 'noAr':
                    src = heatmap(geneList, currentSamples)
                else:
                    src = heatmapAr(geneList, currentSamples, ar)
            context = {'plot': src, 'geneList': json.dumps(geneList), 'searchGeneNames': json.dumps(newGeneNames),
                       'currentSamples':json.dumps(currentSamples), 'addableSamples':json.dumps(rmSampleList),
                       'plotType':plotType}
            return render(request, 'plot.html', context)
        else:
            if len(geneList) != 0 and len(currentSamples) ==0:
                info = "No sample in the list"
            elif len(geneList) == 0 and len(currentSamples) !=0:
                info = "No gene in the list"
            else:
                info = "No gene and sample in the list"
            return render(request, 'plot.html', {'noGene': info, 'geneList': json.dumps(geneList),
                                                 'searchGeneNames': json.dumps(newGeneNames),
                                                 'currentSamples':json.dumps(currentSamples),
                                                 'addableSamples':json.dumps(rmSampleList),
                                                 'plotType':plotType})

def subplot(request):
    if request.method == 'POST':
        #chooseType is the function of adding a gene, removing a gene, adding a sample
        # and removing a sample
        if 'chooseType' in request.POST:
            global ar
            global currentSamples
            ar = request.POST.get('changeGroupBy')
            if len(geneList) != 0 and len(currentSamples) != 0:
                if plotType == 'boxplotT':
                    src = boxplot(geneList, ar, currentSamples)
                elif plotType == 'heatmapT':
                    if ar == 'noAr':
                        src = heatmap(geneList, currentSamples)
                    else:
                        src = heatmapAr(geneList, currentSamples, ar)
                rmSampleList = list(set(currentSamples)^set(allSamples))
                context = {'plot': src, 'geneList': json.dumps(geneList),
                           'searchGeneNames': json.dumps(newGeneNames),
                            'currentSamples':json.dumps(currentSamples),
                           'addableSamples':json.dumps(rmSampleList),
                           'plotType':plotType}
                return render(request, 'plot.html', context)
        #selectSamples is the function of select specific samples or all samples
        elif 'selectSamples' in request.POST:
            if len(geneList) != 0:
                currentSamples = request.POST.getlist('sampleCheckBoxList')
                all = request.POST.getlist('allCheckBoxList')
                #Checking whether user chooses All Samples. If so, current samples will be all samples
                if len(all) != 0:
                    currentSamples = ['1672MH8C','1672MX24C','1672MY25C', '1672MI9T','1672MK11T','1672MM13T', '2011AC3C',
                                      '2011AF6C','2011AG7C','2011AH8C','2011AA1T','2011AB2T', '2011AD4T', '2011AE5T',
                                      '272AG7C','272AP16C','272AW23C','272AZ26C','272AB2T','272ACC29T', '272ADD30T',
                                      '272AS19T', '2012ADD30C','2012AF6C','2012AR18C','2012AX24C','2012AD4T','2012AJ10T',
                                      '2012AW23T','2012AY25T']
                rmSampleList = list(set(currentSamples) ^ set(allSamples))
                if plotType == 'heatmapT':
                    if ar == 'noAr':
                        img = heatmap(geneList, currentSamples)
                    else:
                        img = heatmapAr(geneList, currentSamples, ar)
                elif plotType == 'boxplotT':
                    img = boxplot(geneList, ar, currentSamples)
                return render(request, 'plot.html', context={'plot': img, 'geneList': json.dumps(geneList),
                               'searchGeneNames': json.dumps(newGeneNames), 'addableSamples':json.dumps(rmSampleList),
                               'currentSamples':json.dumps(currentSamples), 'plotType':plotType})
            else:
                info = 'No gene in the list'
                return render(request, 'plot.html', context={'noGene': info})

#Autocomplete searching function.
class GeneNameAutocomplete(autocomplete.Select2QuerySetView):
    def get_queryset(self):
        # Don't forget to filter out results depending on the visitor !
        if not self.request.user.is_authenticated:
            return Geneexpression.objects.none()

        qs = Geneexpression.objects.all()

        if self.q:
            qs = qs.filter(name__istartswith=self.q)
        return qs

#Converting a list to a string. A supporting function for Plot.
def convertListToString(org_list, seperator=','):
    """ Convert list to string, by joining all item in list with given separator.
        Returns the concatenated string """
    return seperator.join(org_list)
