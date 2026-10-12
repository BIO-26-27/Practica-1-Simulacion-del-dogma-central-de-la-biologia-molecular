def transcribir_adn(adn):

    print("\n" + "=" * 60)
    print("TRANSCRIPCIÓN")
    print("=" * 60)

    print("\n1. ARN POLIMERASA")
    print("La ARN polimerasa utiliza la cadena molde de ADN")
    print("para sintetizar una molécula de ARN mensajero.")

    # La secuencia que tenemos corresponde a la cadena codificante.
    cadena_codificante = adn

    # Obtenemos la cadena molde.
    cadena_molde = adn.complement()

    print("\nCadena codificante:")
    print("5'", cadena_codificante[:60], "... 3'")

    print("\nCadena molde:")
    print("3'", cadena_molde[:60], "... 5'")

    # Biopython realiza la transcripción.
    arn_m = cadena_codificante.transcribe()

    print("\n2. ARN MENSAJERO")
    print("Durante la transcripción, la timina (T) es sustituida por uracilo (U).")

    print("\nARNm:")
    print("5'", arn_m[:60], "... 3'")

    print("\nLongitud del ARNm:", len(arn_m), "nucleótidos")

    return arn_m