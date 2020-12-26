import os
import shutil

def write_file_py(dataset, path):
    try:
        filename, file_extension = os.path.splitext(path)
        if (file_extension=='.zip'):
            cmd='mkdir datasets/user_uploaded/'+dataset+'_10x'
            os.system(cmd)
            shutil.move(path, "datasets/user_uploaded/" + dataset+'_10x/'+dataset
                        + '.zip')
            os.chdir('datasets/user_uploaded/'+ dataset+'_10x')
            os.system('unzip '+dataset+'.zip')
            os.chdir('../../..')
        elif (file_extension=='.gz'):
            cmd='mkdir datasets/user_uploaded/'+dataset+'_10x'
            os.system(cmd)
            shutil.move(path, "datasets/user_uploaded/" + dataset+'_10x/'+dataset
                        + '.tar.gz')
            os.chdir('datasets/user_uploaded/'+ dataset+'_10x')
            os.system('tar zxvf '+dataset+'.tar.gz')
            os.chdir('../../..')
        else:
            shutil.move(path, "datasets/user_uploaded/" + dataset + file_extension)
        
    except Exception as e:
        print(str(e))
        raise OSError("A problem occured when reading dataset.")
