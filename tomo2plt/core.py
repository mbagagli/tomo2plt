from pathlib import Path
import numpy as np
import time

from obspy.geodetics import locations2degrees, gps2dist_azimuth
from obspy.geodetics.base import kilometers2degrees, degrees2kilometers
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable


def cpt2cmap(filename, n_bins=256):
    if not filename:
        raise ValueError("I need a valid *cpt file path!")
    #
    with open(filename, 'r') as f:
        lines = f.readlines()

    x = []
    r = []
    g = []
    b = []

    for line in lines:
        if line[0] != '#':
            lf = line.strip().split()
            x.append(float(lf[0]))
            r.append(float(lf[1])/255.)
            g.append(float(lf[2])/255.)
            b.append(float(lf[3])/255.)

    min_x = min(x)
    max_x = max(x)

    x = [(i-min_x)/(max_x-min_x) for i in x]

    cdict = {'red': [], 'green': [], 'blue': []}

    for i in range(len(x)):
        cdict['red'].append([x[i], r[i], r[i]])
        cdict['green'].append([x[i], g[i], g[i]])
        cdict['blue'].append([x[i], b[i], b[i]])

    return mcolors.LinearSegmentedColormap('my_colormap', cdict, n_bins)


def psxy(inax, file_path, delimiter=">", first_is_reference=False):
    """ Read a txt file simil to GMT for line segment """
    # ------ 1. Read file
    def __parse_file__(file_path, delimiter=">"):
        data_dict = {}
        current_section = None
        current_matplotlib_par = {}
        current_plot_array = []

        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()

                if not line or line.startswith('#'):
                    continue  # Ignore empty lines and lines starting with #

                if line.startswith(delimiter):
                    if current_section is not None:
                        # Store data for previous section
                        data_dict[current_section] = {
                            'plt_par': current_matplotlib_par,
                            'plt_array': np.array(current_plot_array)
                        }

                    # Start a new section
                    current_section = f"Sect_{len(data_dict) + 1}"
                    current_matplotlib_par = {}
                    current_plot_array = []

                    # Extract matplotlib parameters
                    matplotlib_params = line[1:].strip().split(',')
                    for param in matplotlib_params:
                        key, value = param.split('=')
                        current_matplotlib_par[key.strip()] = value.strip()
                else:
                    # Parse plot points
                    coords = [float(x) for x in line.split()]
                    current_plot_array.append(coords)

        # Store data for the last section
        data_dict[current_section] = {
            'plt_par': current_matplotlib_par,
            'plt_array': np.array(current_plot_array)
        }
        return data_dict
    #
    plot_dict = __parse_file__(file_path, delimiter=delimiter)
    # breakpoint()
    for kk, vv in plot_dict.items():
        if first_is_reference:
            _refx, _refy = vv['plt_array'][0]
            for row in vv['plt_array']:
                inax.plot(
                    [_refx, row[0]], [_refy, row[1]], **vv['plt_par'], zorder=1)
        else:
            inax.plot(vv['plt_array'][:, 0], vv['plt_array'][:, 1], **vv['plt_par'])
    #
    return inax


class TomoGrid:
    def __init__(self, file_path, **kwargs):
        self.data = self.__load_grid__(file_path, **kwargs)
        self.coordinates = self.data[:, [0, 1, 4]]
        self.xyz = self.data[:, [2, 3, 4]]
        #
        self.pvel = self.data[:, 5]
        self.pdel = self.data[:, 6]
        self.prde = self.data[:, 7]
        self.khit = self.data[:, 8]
        self.dws = self.data[:, 9]
        #
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        string_format = (  # 0          1      2  3    4       5       6
                "Columns: LONGITUDE, LATITUDE, X, Y, DEPTH, Velocity, Delta, "
                # 7    8     9
                "RDE, KHIT, DWS\n"
                "Total Nodes:  %d" % self.data.shape[0]
            )
        return string_format

    def __load_grid__(self, file_path, **kwargs):
        fp = Path(file_path)
        #
        if fp.exists():
            # Change delimiter according to your file
            print("Reading:  %s" % fp.name)
            _data = np.genfromtxt(str(fp.resolve()))  # , delimiter=" ")
        else:
            raise ValueError("FILE:  %s  missing!" % str(fp.resolve()))
        #
        return _data

    def __mode_interpolate__(self, plane, step_x, step_y, what, interpolate):
        _sss = time.time()
        if what.lower() in ("p", "pvel", "p-velocity"):
            values_on_plane = griddata(self.coordinates, self.pvel, plane,
                                       method=interpolate).reshape(
                                            step_x, step_y)

        elif what.lower() in ("delta", "min1d"):
            values_on_plane = griddata(self.coordinates, self.pdel, plane,
                                       method=interpolate).reshape(
                                            step_x, step_y)

        elif what.lower() in ("res", "resolution", "rde"):
            values_on_plane = griddata(self.coordinates, self.prde, plane,
                                       method=interpolate).reshape(
                                            step_x, step_y)

        else:
            raise ValueError("Only 'RES', 'PVEL', or 'DELTA' allowed!")

        _eee = time.time()
        print("Total interpolation TIME:  %.2f min." % ((_eee-_sss)/60.0))
        return values_on_plane

    # ============================================================
    # ============================================================  DEPTH
    # ============================================================
    # ============================================================

    def get_depth_slice(self, start, end,
                        increment_x=0.09, increment_y=0.09,
                        depth=20, what="pvel", interpolate="linear"):

        start = np.array(start)
        end = np.array(end)
        num_steps_x = int(np.abs(start-end)[0]/increment_x)
        num_steps_y = int(np.abs(start-end)[1]/increment_y)

        x = np.linspace(start[0], end[0], num_steps_x)
        y = np.linspace(start[1], end[1], num_steps_y)
        z = np.full((num_steps_x, num_steps_y), depth)
        x_grid, y_grid = np.meshgrid(x, y)
        x_grid = x_grid.T
        y_grid = y_grid.T
        plane_points = np.dstack([x_grid, y_grid, z]).reshape(-1, 3)

        print("Extracting Depth SLICE:")
        print("    LON:  %8.4f %8.4f" % (start[0], end[0]))
        print("    LAT:  %8.4f %8.4f" % (start[1], end[1]))
        print("  SHAPE:  (%d, %d)" % (plane_points.shape[0],
                                      plane_points.shape[1]))

        values_on_plane = self.__mode_interpolate__(
                                    plane_points, num_steps_x, num_steps_y,
                                    what, interpolate)
        print(" ... done!")
        return (x_grid, y_grid, values_on_plane)

    def plot_depth_slice(self, start, end,
                         increment_x_km=10.0, increment_y_km=10.0,
                         depth=20,  what="delta",
                         smooth=0.9, mask_rde=0.2, mask_rde_alpha=1.0,
                         interpolate="linear", palettes=None,
                         isolines=[]):

        """ Extract and plot """
        (x_grid, y_grid, values_on_plane) = self.get_depth_slice(
                        start, end,
                        increment_x=kilometers2degrees(increment_x_km),
                        increment_y=kilometers2degrees(increment_y_km),
                        what=what, depth=depth,
                        interpolate=interpolate)
        if mask_rde:
            (x_grid_rde,
             y_grid_rde,
             values_on_plane_rde) = self.get_depth_slice(
                            start, end,
                            increment_x=kilometers2degrees(5.0),
                            increment_y=kilometers2degrees(5.0),
                            what="rde", depth=depth,
                            interpolate="nearest")

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)

        cmap = cpt2cmap(palettes)
        if smooth:
            from scipy import ndimage
            values_on_plane = ndimage.gaussian_filter(values_on_plane,
                                                      sigma=smooth)
        #
        asd = plt.pcolormesh(x_grid, y_grid, values_on_plane,
                             cmap=cmap, vmin=-10, vmax=10,
                             shading='gouraud', rasterized=True)
        asd.set_linewidth(0)

        cbar = plt.colorbar(asd, label='Vp (%)', ax=ax, shrink=0.4)
        cbar.set_ticks([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])  # Specify the values where you want the ticks
        cbar.set_ticklabels(['-10', '-8', '-6', '-4', '-2', '0',
                             '2', '4', '6', '8', '10'])

        plt.xlabel('longitude')
        plt.ylabel('latitude')

        # ----------------------------  ISOLINES
        if isolines:
            contours = ax.contour(
                        x_grid, y_grid, values_on_plane,
                        levels=isolines, colors="darkgray",  # colors='black',
                        linestyles='-', linewidths=0.5)

            # Add labels to the contours
            contours.clabel(inline=True, fontsize=9, colors='black')

        # ----------------------------  MASK WITH WHITE
        lower_values = values_on_plane_rde < 0.01+0.001
        upper_values = values_on_plane_rde >= 0.01+0.001

        masked_data = np.where(upper_values, values_on_plane_rde, 1)
        masked_data = np.where(lower_values, masked_data, np.nan)

        # # Plot the masked data with transparency
        c_mask = ax.pcolormesh(x_grid_rde,
                               y_grid_rde,
                               masked_data,
                               cmap='Greys',
                               edgecolors='none',
                               shading='auto',
                               rasterized=True)
        c_mask.set_facecolor('white')

        if mask_rde:
            # Create mask for lower and upper values
            lower_values = values_on_plane_rde < mask_rde+0.001
            upper_values = values_on_plane_rde >= mask_rde+0.001

            masked_data = np.where(upper_values, values_on_plane_rde, 1)
            masked_data = np.where(lower_values, masked_data, np.nan)

            # # Plot the masked data with transparency
            c_mask = ax.pcolormesh(x_grid_rde,
                                   y_grid_rde,
                                   masked_data,
                                   # alpha=0.0,  # fully transparent
                                   # alpha=0.5,
                                   alpha=mask_rde_alpha,
                                   cmap='Greys',
                                   edgecolors='none',
                                   shading='auto',  #shading='gouraud',
                                   rasterized=True)

            # Draw contour line at threshold
            contour = ax.contour(x_grid_rde, y_grid_rde, values_on_plane_rde,
                                 levels=[mask_rde],
                                 linewidths=1.2,
                                 colors='black')

        return (fig, ax, cbar)

    # ============================================================
    # ============================================================  VERTICAL
    # ============================================================
    # ============================================================

    def get_depth_section(self, start, end,
                          increment_x_km=10.0, increment_y_km=10.0,
                          what="pvel", interpolate="linear"):

        start, end = np.array(start), np.array(end)
        latlon_start, latlon_end = np.flipud(start[:2]), np.flipud(end[:2])
        dist_km = degrees2kilometers(
                    locations2degrees(*latlon_start, *latlon_end)
                )
        num_steps_xy = int(dist_km / increment_x_km)
        num_steps_z = int(np.abs(start[2]-end[2]) / increment_y_km)

        x = np.linspace(start[0], end[0], num_steps_xy)
        y = np.linspace(start[1], end[1], num_steps_xy)
        z = np.linspace(start[2], end[2], num_steps_z)
        dist_km_points = np.linspace(0, dist_km, num_steps_xy)
        profile_points = np.column_stack([x, y])

        # Repeat the vector X to match the number of rows in M
        z_repeated = np.tile(z, (profile_points.shape[0], 1))

        # Repeat each row of M P times
        profile_points_repeated = np.repeat(profile_points, z.shape[0],
                                            axis=0)

        # Reshape the repeated X vector to match the shape of M_repeated
        z_reshaped = np.reshape(
                        z_repeated,
                        (profile_points.shape[0] * z.shape[0], -1))

        # Concatenate M_repeated and X_reshaped along the second axis
        plane_points = np.concatenate(
                        (profile_points_repeated, z_reshaped), axis=1)

        print("Extracting Depth SECTION @ %s:" % what)
        print("    START:  %8.4f %8.4f" % (start[0], start[1]))
        print("      END:  %8.4f %8.4f" % (end[0], end[1]))
        print("    SHAPE:  (%d, %d)" % (plane_points.shape[0],
                                        plane_points.shape[1]))

        values_on_plane = self.__mode_interpolate__(
                                    plane_points, num_steps_xy, num_steps_z,
                                    what, interpolate)
        x_grid, y_grid = np.meshgrid(dist_km_points, z)
        x_grid = x_grid.T
        y_grid = y_grid.T

        print(" ... done!")
        return (x_grid, y_grid, dist_km, values_on_plane)

    def plot_depth_section(self, start, end,
                           increment_x_km=10.0, increment_y_km=10.0,
                           depth=20,  what="pvel", interpolate="linear",
                           smooth=False, mask_rde=False, mask_rde_alpha=1.0,
                           add_topography=False, palettes=None,
                           isolines=[]):

        """ Extract and plot """
        (x_grid, y_grid, profile_km, values_on_plane) = self.get_depth_section(
                        start, end,
                        increment_x_km=increment_x_km,
                        increment_y_km=increment_y_km,
                        what=what,
                        interpolate=interpolate)
        if mask_rde:
            (x_grid_rde,
             y_grid_rde,
             profile_km,
             values_on_plane_rde) = self.get_depth_section(
                            start, end,
                            increment_x_km=1,
                            increment_y_km=1,
                            what="rde",
                            interpolate="nearest")

        # -------------------------------------------------------------
        # --------------------------------------------  PEREPARE FIGURE

        aspect_ratio = 2
        fixed_height_cm = 7

        # Function to calculate required width based on x-axis limits and aspect ratio
        def calculate_required_width(x_min, x_max):
            return fixed_height_cm * aspect_ratio * (x_max - x_min) / 100

        # Calculate required width for Figure 1
        x_min1, x_max1 = 0, profile_km
        width1 = calculate_required_width(x_min1, x_max1)

        # Convert fixed height and calculated widths to inches
        fixed_height_inches = fixed_height_cm / 2.54
        width1_inches = width1 / 2.54

        fig = plt.figure(figsize=(width1_inches, fixed_height_inches))
        ax1 = fig.add_subplot(1, 1, 1)

        # -------------------------------------------------------------
        # --------------------------------------------  TOPOGRAPHY

        if add_topography:
            DEM = Plane_2D_Grid(add_topography, tag="DEM")
            print("Plotting DEM:  %s" % add_topography)
            (_topo_prof, _topo_height) = DEM.get_profile(
                                                start, end, 2.0,
                                                query_profile_only=True)
            _topo_height[_topo_height < 0.0] = 0.0  # np.nan

            ax1.plot(_topo_prof, -(_topo_height/1000.0)*4, color="black",
                     zorder=9999)
            # To remove the box around the subplot
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

        # -------------------------------------------------------------
        # ---------------------------------------------  VEL.MODEL

        if smooth:
            from scipy import ndimage
            values_on_plane = ndimage.gaussian_filter(values_on_plane,
                                                      sigma=smooth)
        #
        cmap = cpt2cmap(palettes)
        asd = plt.pcolormesh(x_grid, y_grid, values_on_plane,
                             cmap=cmap, vmin=2.5, vmax=9, edgecolors='none',
                             shading='gouraud', rasterized=True)
        asd.set_linewidth(0)

        # ---------------------------- Colorbar
        divider = make_axes_locatable(ax1)
        cbar_ax = divider.append_axes("left", size="1.3%", pad=1.3)
        cbar = plt.colorbar(asd, label='Vp (km/s)',
                            shrink=0.9, cax=cbar_ax)
        cbar.set_ticks([3, 4, 5, 6, 7, 8, 9])  # Specify the values where you want the ticks
        cbar.set_ticklabels(['3', '4', '5', '6', '7', '8', '9'])
        cbar.ax.invert_yaxis()

        # ----------------------------  ISOLINES
        if isolines:
            contours = ax1.contour(
                        x_grid, y_grid, values_on_plane,
                        levels=isolines, colors="darkgray",  # colors='black',
                        linestyles='-', linewidths=0.5)

            # Add labels to the contours
            contours.clabel(inline=True, fontsize=9, colors='darkgray')

        # ----------------------------  7.25 Tomographic MOHO
        ax1.contour(x_grid, y_grid, values_on_plane,
                    levels=[7.25], colors='white',
                    linestyles='dashed', linewidths=1.3)

        # ----------------------------  MASK WITH WHITE
        lower_values = values_on_plane_rde < 0.01+0.001
        upper_values = values_on_plane_rde >= 0.01+0.001

        masked_data = np.where(upper_values, values_on_plane_rde, 1)
        masked_data = np.where(lower_values, masked_data, np.nan)

        # # Plot the masked data with transparency
        c_mask = ax1.pcolormesh(x_grid_rde,
                                y_grid_rde,
                                masked_data,
                                cmap='Greys',
                                edgecolors='none',
                                shading='auto',
                                rasterized=True, zorder=9998)
        c_mask.set_facecolor('white')

        # ----------------------------  Resolution Matrix
        if mask_rde:
            # Create mask for lower and upper values
            lower_values = values_on_plane_rde < mask_rde+0.001
            upper_values = values_on_plane_rde >= mask_rde+0.001
            masked_data = np.where(upper_values, values_on_plane_rde, 1)
            masked_data = np.where(lower_values, masked_data, np.nan)

            # Plot the masked data with transparency
            c_mask = ax1.pcolormesh(x_grid_rde,
                                    y_grid_rde,
                                    masked_data,
                                    alpha=mask_rde_alpha, cmap='Greys',
                                    edgecolors='none', shading='auto',  #shading='gouraud',
                                    rasterized=True)

            # Draw contour line at threshold
            contour = ax1.contour(x_grid_rde, y_grid_rde, values_on_plane_rde,
                                  levels=[mask_rde], colors='black',
                                  linewidths=1.2,
                                  linestyles='dashed')

        # ----------------------------  Finish Axis
        ax1.set_xlabel('distance along profile (km)')
        ax1.set_ylabel('depth (km)')
        ax1.set_xlim([0, np.max(x_grid)])
        ax1.set_ylim([-20, np.max(y_grid)])
        ax1.invert_yaxis()
        ax1.set_yticks(list(np.arange(-20, np.max(y_grid)+5, 10)))
        ax1.set_yticklabels([
                    "-5", "-2.5", "0", "10", "20",
                    "30", "40", "50", "60", "70"])
        ax1.set_aspect(2, adjustable="box", anchor="SW")

        return (fig, ax1, cbar)


class Plane_2D_Grid:
    """ Necessary to extract profiles of value from Start to End.
        Valid for X Y VALUE
    """

    def __init__(self, file_path, tag="xyzGrid", **kwargs):
        if isinstance(file_path, (str, Path)):
            self.data = self.__load_grid__(file_path, **kwargs)
        elif isinstance(file_path, np.ndarray):
            assert file_path.shape[1] == 3
            self.data = file_path
        else:
            raise ValueError("Unsupported input type:  %s" % type(file_path))
        #
        self.coordinates = self.data[:, :-1]
        self.values = self.data[:, -1]
        self.tag = tag
        #
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __load_grid__(self, file_path, **kwargs):
        fp = Path(file_path)
        if fp.exists():
            _data = np.loadtxt(str(fp.resolve()), **kwargs)  #, delimiter=" ")  # Change delimiter according to your file
        else:
            raise ValueError("FILE:  %s  missing!" % str(fp.resolve()))
        return _data

    def get_profile(self, start, end, step_km, query_profile_only=False):
        """ Extract profile and returns the """
        start, end = np.array(start), np.array(end)
        latlon_start, latlon_end = np.flipud(start[:2]), np.flipud(end[:2])
        dist_km = degrees2kilometers(
                    locations2degrees(*latlon_start, *latlon_end)
                )

        num_steps_xy = int(dist_km / step_km)

        if query_profile_only:
            # Extract only important parts from the grid
            # 0. define bound
            lower_lon_bound = min(start[0], end[0])
            upper_lon_bound = max(start[0], end[0])
            lower_lat_bound = min(start[1], end[1])
            upper_lat_bound = max(start[1], end[1])

            if lower_lon_bound == upper_lon_bound:
                _plot_data = self.data[
                                    (self.data[:,0] >= lower_lon_bound-0.1) &
                                    (self.data[:,0] <= upper_lon_bound+0.1) &
                                    (self.data[:,1] >= lower_lat_bound) &
                                    (self.data[:,1] <= upper_lat_bound)]
                _plot_coord, _plot_values = _plot_data[:,:2], _plot_data[:, 2]

            elif lower_lat_bound == upper_lat_bound:
                _plot_data = self.data[
                                    (self.data[:,0] >= lower_lon_bound) &
                                    (self.data[:,0] <= upper_lon_bound) &
                                    (self.data[:,1] >= lower_lat_bound-0.1) &
                                    (self.data[:,1] <= upper_lat_bound+0.1)]
                _plot_coord, _plot_values = _plot_data[:,:2], _plot_data[:, 2]

            else:
                _plot_data = self.data[
                                    (self.data[:,0] >= lower_lon_bound) &
                                    (self.data[:,0] <= upper_lon_bound) &
                                    (self.data[:,1] >= lower_lat_bound) &
                                    (self.data[:,1] <= upper_lat_bound)]
                _plot_coord, _plot_values = _plot_data[:,:2], _plot_data[:, 2]

        else:
            _plot_coord, _plot_values = self.coordinates, self.values

        if _plot_coord.size == 0 and _plot_values.size == 0:
            raise ValueError("NO POINTS to query on !!!")

        x = np.linspace(start[0], end[0], num_steps_xy)
        y = np.linspace(start[1], end[1], num_steps_xy)

        dist_km_points = np.linspace(0, dist_km, num_steps_xy)
        profile_points = np.column_stack([x, y])

        # Interpolate
        _sss = time.time()
        profile_values = griddata(_plot_coord,
                                  _plot_values,
                                  profile_points, method="linear")
        _eee = time.time()
        print("Total interpolation TIME:  %.2f min" % ((_eee-_sss)/60.0))

        return (dist_km_points, profile_values)


class Points_2D_Grid:
    """ Necessary to extract projected points
        Valid for X Y VALUE
    """

    def __init__(self, file_path, tag="xyzGrid", **kwargs):
        if isinstance(file_path, (str, Path)):
            self.data = self.__load_grid__(file_path, **kwargs)
        elif isinstance(file_path, np.ndarray):
            assert file_path.shape[1] == 3
            self.data = file_path
        else:
            raise ValueError("Unsupported input type:  %s" % type(file_path))
        #
        self.coordinates = self.data[:, :-1]
        self.values = self.data[:, -1]
        self.tag = tag
        #
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __load_grid__(self, file_path, **kwargs):
        fp = Path(file_path)
        if fp.exists():
            _data = np.loadtxt(str(fp.resolve()), **kwargs)  #, delimiter=" ")  # Change delimiter according to your file
        else:
            raise ValueError("FILE:  %s  missing!" % str(fp.resolve()))
        return _data

    def project_points(self, start, end, project=None):
        """ project args is in km """
        # -----------------------------  Inner methods
        def calculate_distance_to_line(point, line_start, line_end):
            # Calculate the distances from the point to the line start and end
            dist_start = locations2degrees(line_start[1],  # lat1
                                           line_start[0],  # lon1
                                           point[1],  # lat2
                                           point[0])  # lon2
            dist_end = locations2degrees(line_end[1],  # lat1
                                         line_end[0],  # lon1
                                         point[1],  # lat2
                                         point[0])  # lon2

            # Calculate the azimuths
            _, az, _ = gps2dist_azimuth(line_start[1],  # lat1
                                        line_start[0],  # lon1
                                        point[1],  # lat2
                                        point[0])  # lon2
            _, az_line, _ = gps2dist_azimuth(line_start[1],  # lat1
                                             line_start[0],  # lon1
                                             line_end[1],  # lat2
                                             line_end[0])  # lon2

            # Calculate perpendicular distance to the line
            dist_to_line = np.sin(np.abs(az - az_line) * np.pi / 180) * min(dist_start, dist_end)
            return degrees2kilometers(dist_to_line)

        def calculate_distance_along_line(point, line_start):
            return degrees2kilometers(locations2degrees(line_start[1],  # lat1
                                                        line_start[0],  # lon1
                                                        point[1],  # lat2
                                                        point[0]))  # lon2
        # --------------------------------------------

        # Set up
        start, end = np.array(start), np.array(end)
        if not project:
            project = 10000  # almost infinite values to project events

        # Calculate and project
        projected_events = []
        for event in self.data:
            dist_to_profile = calculate_distance_to_line(event,
                                                         start[:2],
                                                         end[:2])

            if dist_to_profile <= project:
                distance_along_profile = calculate_distance_along_line(event,
                                                                       start)
                projected_events.append((distance_along_profile, event[2]))  # distance along profile and depth

        projected_events = np.array(projected_events)
        return (projected_events[:, 0], projected_events[:, 1])
