#!/bin/bash

#executar como ./likwid <core> <grupo>
#executa o programa matriz_desempenho com matrizes de tamanho "inicio" "fim"
#Deve-se passar o core a ser utilizado e indicador na linha de comando (L2CACHE, L3, FLOPS_DP e FLOPS_AVX)
#o resultado é um arquivo com 5 colunas: ordem da matriz, resultados do indicador para cada uma 
#das 4 funções

LIKWID_CMD="likwid-perfctr -C $1 -g $2 -m -O"

cache="L2 miss ratio"
cachebw="L3 bandwidth"
flopsDP="Scalar MUOPS/s"
flopsAVX="packed MUOPS/s"

ordemMatriz="$3"

#tamanho inicial e final da matriz, sendo que de inicio a final vai incrementando em base 2
inicio=8
fim=2048
#passo e o multiplicador para os numeros ficarem na base 2
passo=2

#fazer isso se for executar no dinf :)
#export PATH=/home/soft/likwid/bin:/home/soft/likwid/sbin:$PATH
#export LD_LIBRARY_PATH=/home/soft/likwid/lib:$LD_LIBRARY_PATH

#fazer isso no pc de casa :)
#modprobe msr

# Para obter lista de grupos de indicadores de performance:
#      likwid-perfctr -a

# Para obter topologia dos processadores
#      likwid-topology -c -g

#---Confere os indicadores e salva o padrao que procuramos em cada um---
if [ "$2" == "L2CACHE" ]
then
	padrao="$cache"
elif [ "$2" == "L3" ]
then
	padrao="$cachebw"
elif [ "$2" == "FLOPS_DP" ]
then
	padrao="$flopsDP"
else
	padrao="$flopsAVX"
fi

#-----------------------------------------------------------------------

while [ "$inicio" -le "$fim" ]
do
	${LIKWID_CMD} ./matriz_desempenho -n $inicio > temp.tmp

	#ordem da matriz em Bytes
	printf "$((inicio*8)) "
	#procuro pelo "padrao" no arquivo de saida da execucao do programa
	#usando o separador "," salvo a segunda coluna ($2)
	#troco "\n" por " "
	printf "$(grep "$padrao" temp.tmp | awk -F"," '{print $2}' | tr "\n" " ")\n"

	(( inicio = inicio*"$passo" ))

#ao final do while, salvo tudo em um arquivo com o nome do padrao
done > $2.tmp

#-----------------------------------------------------------------------
#----------------gnuplot para fazer os graficos-------------------------
gnuplot <<- EOF
	set logscale x
	set xlabel "Ordem Matriz (bytes)"
	set ylabel "$escolha"
	set title "Medição de performance para $2"   
	set style data point
	set style function line

	set style line 1 lc rgb "red" lw 2
	set style line 2 lc rgb "orange" lw 2
	set style line 3 lc rgb "green" lw 2
	set style line 4 lc rgb "blue" lw 2

	set term png
	set output "$2.png"
	plot "$2.tmp" using 1:2 ls 1 title 'multMatrizNormal' with lines, \
	"$2.tmp" using 1:3 ls 2 title 'multMatrizTransposta' with lines, \
	"$2.tmp" using 1:4 ls 3 title 'multMatrizNormalBloco' with lines, \
	"$2.tmp" using 1:5 ls 4 title 'multMatrizTranspostaBloco' with lines
EOF

rm *.tmp