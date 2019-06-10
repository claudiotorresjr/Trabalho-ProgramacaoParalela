#!/bin/bash

#pograma a ser rodado recebe o primeiro argumento do script
#programa="$1"

#tamanho inicial e final da matriz, sendo que de inicio a final vai incrementando em base 2
#inicio=8
#fim=4096
#passo e o multiplicador para os numeros ficarem na base 2
#passo=2
#resultados="resultadosBase"
#rm -r "$resultados"

#se nao existe o diretorio "resultadosBase", crio um e calculo varios resultados
#de diversos sistemas lineares de "inicio" ate "fim" tamanho
if [[ ! -e "$resultados" ]] 
then
	mkdir "$resultados"

	for proc in 1 #2 4 8 16 32
	do
		for tam in 3000 4243 6000 8486 12000
		do
			for exec in {1..20}
			do
				sudo nice -n -20 ./LU_base/sisLinear -n "$tam" > "$resultados"/"$tam".temp
				soma=$(grep "Time" "$resultados"/"$tam".temp | awk -F":" '{print $2}')
				echo "$soma" >> soma.tmp
			done
			#numero de processadores
			printf "$proc "
			#quantidade de elementos da matriz em mbytes
			printf "$((tam*tam*8/1024/1024)) "
			#printf "$((tam)) "
			#media e desvio padrao do tempo de execucao
			./desvio 
			rm *.tmp
		done 
	done > "$resultados"/Base.txt
fi

#-----------------------------------------------------------------------
#----------------gnuplot para fazer os graficos-------------------------
gnuplot <<- EOF
	set logscale x
	set xlabel "Ordem Matriz (Mbytes)"
	set ylabel "Tempo em segundos"
	set title "Decomposicao LU Base"   
	set style data point
	set style function line

	set style line 1 lc rgb "red" lw 2
	set style line 2 lc rgb "orange" lw 2
	set style line 3 lc rgb "green" lw 2
	set style line 4 lc rgb "blue" lw 2

	set term png
	set output "LU_base/LU_base.png"
	plot ""$resultados"/Base.txt" using 1:2 ls 1 title 'multMatrizNormal' with lines, \
	""$resultados"/Base.txt" using 1:3 ls 2 title 'multMatrizTransposta' with lines, \
	""$resultados"/Base.txt" using 1:4 ls 3 title 'multMatrizNormalBloco' with lines, \
	""$resultados"/Base.txt" using 1:5 ls 4 title 'multMatrizTranspostaBloco' with lines
EOF


#versao otimizada
resultados="resultadosLU_otimizado"
if [[ ! -e "$resultados" ]] 
then
	mkdir "$resultados"
	for proc in 1 2 4 8 16 32
	do
		for tam in 3000 4243 6000 8486 12000
		do
			for exec in {1..20}
			do
				sudo nice -n -20 ./LU_otimizado/sisLinear -p "$proc" -n "$tam" > "$tam".temp
				soma=$(grep "Time" "$tam".temp | awk -F":" '{print $2}')
				echo "$soma" >> soma.tmp
			done
			
			#numero de processadores
			printf "$proc " >> "$resultados"/LU_otimizado.txt
			#quantidade de elementos da matriz em bytes
			printf "$((tam*tam*8/1024/1024)) " >> "$resultados"/LU_otimizado.txt
			#media e desvio padrao do tempo de execucao
			./desvio >> "$resultados"/LU_otimizado.txt
			rm *.tmp
		done 
	done
fi
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

rm *.temp


resultados="resultadosLU_paralelo"
if [[ ! -e "$resultados" ]] 
then
	mkdir "$resultados"
	for proc in 1 2 4 8 16 32
	do
		for tam in 3000 4243 6000 8486 12000
		do
			for exec in {1..20}
			do
				sudo nice -n -20 ./LU_paralelo/sisLinear -p "$proc" -n "$tam" > "$tam".temp
				soma=$(grep "Time" "$tam".temp | awk -F":" '{print $2}')
				echo "$soma" >> soma.tmp
			done
			
			#numero de processadores
			printf "$proc " >> "$resultados"/LU_paralelo.txt
			#quantidade de elementos da matriz em bytes
			printf "$((tam*tam*8/1024/1024)) " >> "$resultados"/LU_paralelo.txt
			#media e desvio padrao do tempo de execucao
			./desvio >> "$resultados"/LU_paralelo.txt
			rm *.tmp
		done 
	done
fi
rm *.temp