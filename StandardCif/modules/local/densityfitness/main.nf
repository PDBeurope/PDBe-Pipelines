process DENSITY_FITNESS {
    container 'file:///Users/ahsan/OrbStack/docker/images/density-fitness'

    input:
    path cif_file
    path mtz_file

    output:
    path "myfile*" // Adjust based on the actual extension the tool produces

    script:
    """
    density-fitness \\
        --hklin $mtz_file \\
        --xyzin $cif_file \\
        --output myfile
    """
}