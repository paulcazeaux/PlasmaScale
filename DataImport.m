function [density, velocity, pressure] = DataImport(filename)
    %%% We are going to parse the file 'filename' in order to extract the
    %%% time series for density, velocity and pressure obtained and 
    %%% exported by the PlasmaScale code.
    
    density = timeseries('Ion Density');
    velocity = timeseries('Ion mean velocity');
    pressure = timeseries('Ion Pressure');
    
    fid = fopen(filename, 'r');
    if (fid == -1)
        print('Data file not found');
    end
    
    [t, count] = fscanf(fid, '\nt = %f');
    while (count)
        fscanf(fid, 'Density:');
        dens = fscanf(fid, '%f');
        density = addsample(density, 'Data', dens', 'Time', t);
        fscanf(fid, 'Velocity:');
        vel = fscanf(fid, '%f');
        velocity = addsample(velocity, 'Data', vel', 'Time', t);
        fscanf(fid, 'Pressure:');
        p = fscanf(fid, '%f');
        pressure = addsample(pressure, 'Data', p', 'Time', t);
        [t, count] = fscanf(fid, '\nt = %f');
    end
    fclose(fid);
end