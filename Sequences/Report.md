# Reporte de construcción de biblioteca de secuencias de AVCR

_Erwin Quiroga_

## Introducción

La busqueda de las secuencias se realizó usando BLAST de manera remota, usando el modulo `Bio.BLAST.NCBIWWW`. Los targets fueron obtenidos desde el Uniprot y se realizó el blastp contra la base de datos Swissprot.


|Proteina|Codigo Uniprot|
|  :----:  |  :-----:   |
|  ALK1  |  P37023      |
|ALK2    |Q04771        |
|ALK3    |P36894  |
|ALK4    |P36896|
|ALK5    |P36897|
|ALK6    |Q05438|
|ALK7    |Q8NER5|
|BMPR2   |Q13873|
|ActR-IIA|Q7SXW6|
|ActR-IIB|Q56A35|
|TGFR2   |P37173|
|AMH-RII |Q16671|
|EGFR (control)|P00533|


El Blast se realizo con los siguientes parametros

|  Parametro | Valor  |
|:-:|:-:|
| Programa  |  Blastp |
| Base de datos  |  SwissProt |
| Tolerancia (Corte del E-Value)  | 0.01  |
|  Hitlist | 2000  |

Los resultados de los BLAST fueron graficados en función del logaritmo del E-Value contra su posición en el listado de los hits, siguiendo el siguiente código:


>```python
>
>result_handle = open(blastfile)
>record = NCBIXML.read(result_handle)
>
>evalue=[]
>for aln in record.alignments:
>    hsp = aln.hsps[0]
>    evalue.append(hsp.expect)
>
>log_eval=[]
>for i in evalue:
>    if i == 0:
>        log_eval.append(-400)
>    else:
>        log_eval.append(np.log(i))
>
>d_eval= list(np.diff(log_eval))
>d_eval.insert(0,d_eval[0])
>index = list(range(len(evalue)))
>```

Lo que generó los siguentes graficos, generados con bohek



|![alt text](http://localhost:8888/files/Documentos/AVCR/Sequences/plot_ALK1.png "Blast de ALK1 conta SwissProt")|
|:---:|
|_Blast de ALK1 conta SwissProt_|


