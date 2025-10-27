import os
import pandas as pd

directory = '/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files/primary-exoplanet'

names = [os.path.splitext(file)[0].replace('_target_definition', '') 
         for file in os.listdir(directory) 
         if file.endswith('_target_definition.json')]

# names = ['WASP-107b', 'TOI-1685b', 'WASP-52b', 'TOI-1416b', 'HIP_65_Ab', 'K2-198b', 'HD_3167b', 'WASP-80b', 'TOI-3884b', 'LTT_1445_Ac', 'WASP-177b', 'WASP-69b', 'TOI-942b', 'GJ_9827b', 'GJ_1214b', 'TOI-836b', 'TOI-244b', 'TOI-776b', 'TOI-2427b', 'L_98-59d']

# Create a DataFrame
def format_names(name):
    # Special case for HIP_65_Ab
    if name == 'HIP_65_Ab':
        return 'HIP 65 A b', 'HIP 65 A'
    
    # Special case for LTT_1445_Ac
    if name == 'LTT_1445_Ac':
        return 'LTT 1445 A c', 'LTT 1445 A'
    
    # Special case for L_98-59d
    if name == 'L_98-59d':
        return 'L 98-59 d', 'L 98-59'

    if name == 'HD_3167b':
        return 'HD 3167 b', 'HD 3167'

    if name == 'GJ_9827b':
        return 'GJ 9827 b', 'GJ 9827'

    if name == 'GJ_1214b':
        return 'GJ 1214 b', 'GJ 1214'
    
    # For other names, insert a space before the last character
    planet_name = name[:-1] + ' ' + name[-1]
    
    # Star name is everything before the last character
    star_name = name[:-1].replace('_', ' ')
    
    return planet_name, star_name

# Create lists for DataFrame
planet_names = []
star_names = []

for name in names:
    planet_name, star_name = format_names(name)
    planet_names.append(planet_name)
    star_names.append(star_name)

# Create a DataFrame
df = pd.DataFrame({
    'Planet Name': planet_names,
    'Planet Simbad Name': planet_names,
    'Star Name': star_names,
    'Star Simbad Name': star_names,
    'Number of Transits to Capture': [10 for _ in names],  # Placeholder
    'Priority': [1 for _ in names]  # Placeholder
})

# Reorder columns to match the desired output
df = df[['Planet Name', 'Planet Simbad Name', 'Star Name', 'Star Simbad Name', 'Number of Transits to Capture', 'Priority']]