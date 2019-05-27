#!/bin/bash

#pograma a ser rodado recebe o primeiro argumento do script
programa="$1"

#tamanho inicial e final da matriz, sendo que de inicio a final vai incrementando em base 2
inicio=8
fim=4096
#passo e o multiplicador para os numeros ficarem na base 2
passo=2

resultados="resultadosBase"
#se nao existe o diretorio "resultadosBase", crio um e calculo varios resultados
#de diversos sistemas lineares de "inicio" ate "fim" tamanho
if [[ ! -e "$resultados" ]] 
then
    mkdir "$resultados"

    while [ "$inicio" -le "$fim" ]
	do
		./LU_base/sisLinear -n "$inicio" > "$resultados"/"$inicio".txt

		(( inicio = inicio*"$passo" ))
	done
fi

#executa outra versao do programa e verifica se o resultado eh igual ao da versao base
#while [ "$inicio" -le "$fim" ]
#do
#	./"$programa" -n "$inicio" > "$inicio".temp
#
#	if diff "$resultados"/"$inicio".txt "$inicio".temp
#	then
#		echo "certo"
#	else
#		echo "Deu ruim"
#    	exit 1
#    fi
#
#	(( inicio = inicio*"$passo" ))
#done
#
#rm *.temp
