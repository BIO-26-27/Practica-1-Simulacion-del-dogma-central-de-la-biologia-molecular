from Bio import SeqIO
from replicacion import replicar_adn
from transcripcion import transcribir_adn
from traduccion import traducir_arn

# Leer el archivo GenBank
registro = SeqIO.read(
    "simulador-dogma-central/data/lacZ_Ecoli_NC_000913.3.gb",
    "genbank"
)

# Obtener la secuencia
adn = registro.seq


print("=" * 60)
print("SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR")
print("=" * 60)

print("\nDATOS DE LA SECUENCIA")

print("ID:", registro.id)
print("Descripción:", registro.description)
print("Longitud:", len(adn), "nucleótidos")

print("\nPrimeros 60 nucleótidos:")
print(adn[:60])

print("\nÚltimos 60 nucleótidos:")
print(adn[-60:])

cadena_codificante = adn
cadena_molde = adn.complement()

print("\n" + "=" * 60)
print("DOBLE CADENA DE ADN")
print("=" * 60)

print("\nCadena codificante:")
print("5'", cadena_codificante[:60], "3'")

print("  ", "|" * 60)

print("Cadena molde:")
print("3'", cadena_molde[:60], "5'")

# REPLICACIÓN
cadena_lider, cadena_rezagada = replicar_adn(adn)
# TRANSCRIPCIÓN
arn_m = transcribir_adn(adn)
# TRADUCCIÓN
proteina = traducir_arn(arn_m)