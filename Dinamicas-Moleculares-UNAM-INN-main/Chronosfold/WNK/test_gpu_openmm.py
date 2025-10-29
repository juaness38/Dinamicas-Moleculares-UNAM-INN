import openmm
from openmm import Platform

print("="*60)
print("OPENMM GPU TEST")
print("="*60)

n_platforms = Platform.getNumPlatforms()
print(f"\nPlatforms disponibles: {n_platforms}")

for i in range(n_platforms):
    platform = Platform.getPlatform(i)
    print(f"  {i}: {platform.getName()}")

# Test CUDA
try:
    cuda = Platform.getPlatformByName('CUDA')
    print(f"\n✓ CUDA Platform disponible")
    print(f"  Devices: {cuda.getPropertyDefaultValue('DeviceIndex')}")
except:
    print("\n✗ CUDA Platform NO disponible")

print("="*60)
