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
#if [[ ! -e "$resultados" ]] 
#then
#	mkdir "$resultados"
#
#	for proc in 1 #2 4 8 16 32
#	do
#		for tam in 3000 4243 6000 8486 12000
#		do
#			#for exec in {1..20}
#			#do
#				sudo nice -n -20 ./LU_base/sisLinear -n "$tam" > "$resultados"/"$tam".temp
#				#soma=$(grep "Time" "$resultados"/"$tam".temp | awk -F":" '{print $2}')
#				#echo "$soma" >> soma.tmp
#			#done
#			#numero de processadores
#			#printf "$proc "
#			#quantidade de elementos da matriz em mbytes
#			#printf "$((tam*tam*8/1024/1024)) "
#			#printf "$((tam)) "
#			#media e desvio padrao do tempo de execucao
#			#./desvio 
#			#rm *.tmp
#		done 
#	done #> "$resultados"/Base.txt
#fi

#versao otimizada
#resultados="resultadosLU_otimizado"
#if [[ ! -e "$resultados" ]] 
#then
#	mkdir "$resultados"
#	for proc in 1 2 4 8 16 32
#	do
#		for tam in 3000 4243 6000 8486 12000
#		do
#			for exec in {1..20}
#			do
#				sudo nice -n -20 ./LU_otimizado/sisLinear -p "$proc" -n "$tam" > "$tam".temp
#				soma=$(grep "Time" "$tam".temp | awk -F":" '{print $2}')
#				echo "$soma" >> soma.tmp
#			done
#			
#			#numero de processadores
#			printf "$proc " >> "$resultados"/LU_otimizado.txt
#			#quantidade de elementos da matriz em bytes
#			printf "$((tam*tam*8/1024/1024)) " >> "$resultados"/LU_otimizado.txt
#			#media e desvio padrao do tempo de execucao
#			./desvio >> "$resultados"/LU_otimizado.txt
#			rm *.tmp
#		done 
#	done
#fi
#rm *.temp