function content = readSol(nDim,nEqu,inFile,outFile)
%readSol Read a Kardos solution or mesh file.
%
%Syntax:
%
%            readSol(nDim,nEqu,inFile,outFile)
%  content = readSol(nDim,nEqu,inFile)
%  content = readSol(nDim,nEqu,inFile,outFile)
%
%Arguments:
%
%  content  A structure whose field names and data correspond to the key
%           words and respective data in the Kardos file.
%
%  nDim     Number of spatial dimensions of the problem.
%
%  nEqu     Number of governing equations of the problem.
%
%  inFile   The name or path of the Kardos file to be read.
%
%  outFile  The name or path of the Matlab file (*.mat) to be written. The
%           resulting Matlab file has each field of 'content' as individual
%           data.
%
%Examples:
%
%  % The following blocks (A) and (B) are equivalent, as they both provide
%  % a structure 'S' in the workspace and a file 'solution.mat' in the 
%  % working directory.
%
%  % (A)
%  readSol(2,3,'2DFlow0123.sol','solution.mat');
%  S = load('solution.mat');
%
%  % (B)
%  S = readSol(2,3,'2DFlow0123.sol');
%  save('solution.mat','-struct','S');
%
%Remarks:
%
%  This function can easily be modified to provide for further Kardos 
%  keywords, by adding keyword/format pairs or switch cases. Make sure, 
%  however, that you do not introduce duplicates!
%
%See also: writeSol, save, load
%           
%Sebastian Ullmann
%Numerical Analysis and Scientific Computing
%Technische Universität Darmstadt
%2010/11/24

  if nargout == 0 && nargin < 4
    warning(['Function ''readSol'' has no effect, because no file ' ...
             'is written and no output is used.']);             %#ok<WNTAG>
    return
  end
  
  assert(isnumeric(nDim)&&isscalar(nDim),'Check 1st input argument.');
  assert(isnumeric(nEqu)&&isscalar(nEqu),'Check 2nd input argument.');
  assert(ischar(inFile)&&~iscell(inFile),'Check 3rd input argument.');
  
  % Add custom keywords and format specifiers, here:
  
  f.points               = ['%d:' repmat('%f,',1,nDim-1) '%f'];
  f.circle_midpoints     = ['%d:' repmat('%f,',1,nDim-1) '%f'];
  f.boundary_points      = '%d:%d'; 
  f.point_class          = '%d:%d';  
  f.edge_class           = '(%d,%d):%d';  
  f.circle_edges         = '(%d,%d):%d';  
  f.triangles            = '%d:%d,%d,%d';  
  f.tetrahedra           = '%d:%d,%d,%d,%d';  
  f.boundary_edges       = '(%d,%d):%d';  
  f.boundary_triangles   = '(%d,%d,%d),%d';  
  f.tetrahedron_class    = '%d:%d';  
  f.solution_at_point    = ['%d:' repmat('%f,',1,nEqu-1) '%f'];
  f.ut_solution_at_point = ['%d:' repmat('%f,',1,nEqu-1) '%f'];

  content = struct();
  
  fid = fopen(inFile,'r');
  assert(fid ~= -1,['File ''' inFile ''' not found.']);
  
  while true
    name = fgetl(fid);
    if ~ischar(name) && name==-1
      break
    end
    if isempty(name)
      warning(['Unexpected empty line in Kardos file ''' ...
               inFile '''.']);                                  %#ok<WNTAG>
      continue
    end
    switch name
      case 'name_of_triangulation'
        content.(name) = fgetl(fid);
        assert(isempty(fgetl(fid)),'This line should be empty.');
      case {'no_of_points', ...
            'no_of_boundary_points', ...
            'no_of_edges', ...
            'no_of_boundary_edges', ...
            'no_of_triangles', ...
            'no_of_boundary_triangles', ...
            'no_of_tetrahedra'}
        content.(name) = str2double(fgetl(fid));
        assert(isempty(fgetl(fid)),'This line should be empty.');
      case 'boundary_types'
        i = 1;
        while true
          s = fgetl(fid);
          if isempty(s)
            break;
          end
          content.(name){i} = s;
          nBoundaryEquations = length(s)-strfind(s,':');
          assert(nBoundaryEquations == nEqu,['Number of equations is ' ...
            num2str(nEqu) ', but the statement ' s ' in file '...
            inFile ' implies that the number of boundary ' ...
            'equations is ' num2str(nBoundaryEquations) '.']);
          i = i+1;
        end
      case fields(f)
        nColumns = length(strfind(f.(name),'%'));
        content.(name) = fscanf(fid, [f.(name) '\n'], [nColumns inf]);
      case 'end'
        break;
      otherwise
        warning(['Unexpected key word ''' name '''. ' ...
                 'Check Kardos file ''' inFile ''' and ' ...
                 'm-file ''' mfilename('fullpath') '''. ' ...
                 'Respective data is discarded.']);           %#ok<WNTAG>
        while true
          s = fgetl(fid);
          if isempty(s) || ~ischar(s)
            break;
          end
        end
        if ~ischar(name) && name==-1
          break
        end
    end

  end
  
  fclose(fid);
  
  if ~isempty(who('outFile'))
    save(outFile,'-struct','content');
  end
  
end