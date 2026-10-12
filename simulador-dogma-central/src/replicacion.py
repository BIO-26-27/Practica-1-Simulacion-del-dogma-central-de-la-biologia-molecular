from Bio.Seq import Seq


def replicar_adn(adn, tam_okazaki=500, tam_cebador=10):

    print("\n" + "=" * 60)
    print("REPLICACIÓN DEL ADN")
    print("=" * 60)

    # ----------------------------------------------------
    # 1. HELICASA
    # ----------------------------------------------------
    print("\n1. HELICASA")
    print("La helicasa abre la doble hélice separando las dos cadenas.")

    cadena_superior = adn
    cadena_inferior = adn.complement()

    print("\nCadena superior:")
    print("5'", cadena_superior[:60], "3'")

    print("\nCadena inferior:")
    print("3'", cadena_inferior[:60], "5'")


    # ----------------------------------------------------
    # 2. CADENA LÍDER
    # ----------------------------------------------------
    print("\n2. CADENA LÍDER")
    print("La primasa coloca un cebador de ARN.")
    print("La ADN polimerasa sintetiza la nueva cadena de forma continua.")

    cadena_lider = Seq(str(adn))

    # Cebador de ARN
    cebador_lider = cadena_lider[:tam_cebador].transcribe()

    print("\nCebador ARN:")
    print("5'", cebador_lider, "3'")

    print("\nNueva cadena líder:")
    print("5'", cadena_lider[:60], "... 3'")


    # ----------------------------------------------------
    # 3. CADENA REZAGADA
    # ----------------------------------------------------
    print("\n3. CADENA REZAGADA")
    print("La ADN polimerasa solo puede sintetizar en dirección 5' → 3'.")
    print("Por ello, esta cadena se sintetiza mediante fragmentos de Okazaki.")

    fragmentos_okazaki = []

    # Dividimos la cadena superior en fragmentos
    fragmentos_molde = []

    for i in range(0, len(adn), tam_okazaki):
        fragmento = adn[i:i + tam_okazaki]
        fragmentos_molde.append(fragmento)

    # Los fragmentos se sintetizan desde el extremo opuesto
    for fragmento in reversed(fragmentos_molde):

        nuevo_fragmento = fragmento.reverse_complement()

        fragmentos_okazaki.append(nuevo_fragmento)


    # Mostrar los fragmentos
    for numero, fragmento in enumerate(fragmentos_okazaki, start=1):

        cebador = fragmento[:tam_cebador].transcribe()

        print(f"\nFragmento de Okazaki {numero}")
        print("Cebador ARN:", cebador)
        print("ADN:", fragmento[:50], "...")


    # ----------------------------------------------------
    # 4. ADN LIGASA
    # ----------------------------------------------------
    print("\n4. ADN LIGASA")
    print("La ADN ligasa une los fragmentos de Okazaki.")

    cadena_rezagada = Seq("".join(str(f) for f in fragmentos_okazaki))

    print("\nNueva cadena rezagada:")
    print("5'", cadena_rezagada[:60], "... 3'")


    # ----------------------------------------------------
    # RESULTADO
    # ----------------------------------------------------
    print("\n" + "=" * 60)
    print("RESULTADO DE LA REPLICACIÓN")
    print("=" * 60)

    print("\nLa molécula de ADN original ha dado lugar")
    print("a dos moléculas de ADN mediante replicación semiconservativa.")

    return cadena_lider, cadena_rezagada