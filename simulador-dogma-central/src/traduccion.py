from Bio.Seq import Seq


def traducir_arn(arn_m):

    print("\n" + "=" * 60)
    print("TRADUCCIÓN Y SÍNTESIS DE PROTEÍNAS")
    print("=" * 60)

    print("\n1. INICIO DE LA TRADUCCIÓN")
    print("El ribosoma lee el ARNm en grupos de tres nucleótidos")
    print("denominados codones.")

    print("\nPrimeros codones del ARNm:")

    # Mostrar los primeros 10 codones
    for i in range(0, 30, 3):
        codon = arn_m[i:i + 3]
        aminoacido = codon.translate(table=11)

        print(f"{codon} -> {aminoacido}")

    # ----------------------------------------------------
    # TRADUCCIÓN COMPLETA
    # ----------------------------------------------------

    print("\n2. SÍNTESIS DE LA PROTEÍNA")
    print("Cada codón determina un aminoácido.")

    proteina = arn_m.translate(
        table=11,
        to_stop=True
    )

    print("\nPrimeros 60 aminoácidos:")
    print(proteina[:60])

    print("\nÚltimos 60 aminoácidos:")
    print(proteina[-60:])

    print("\nLongitud de la proteína:")
    print(len(proteina), "aminoácidos")

    # ----------------------------------------------------
    # CODÓN DE TERMINACIÓN
    # ----------------------------------------------------

    codon_stop = arn_m[-3:]

    print("\n3. TERMINACIÓN")
    print("Último codón del ARNm:", codon_stop)

    if str(codon_stop) in ["UAA", "UAG", "UGA"]:
        print("Es un codón STOP.")
        print("El ribosoma finaliza la traducción.")
    else:
        print("No se ha encontrado un codón STOP al final.")

    print("\n" + "=" * 60)
    print("RESULTADO DE LA TRADUCCIÓN")
    print("=" * 60)

    print("\nSe ha sintetizado la proteína a partir del ARNm.")

    return proteina