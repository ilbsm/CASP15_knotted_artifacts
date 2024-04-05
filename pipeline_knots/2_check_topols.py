from topoly import alexander
from tqdm import tqdm
import os
import signal

class timeout:                                                                  
    def __init__(self, seconds=1, error_message='Timeout'):                     
        self.seconds = seconds                                                  
        self.error_message = error_message                                      
    def handle_timeout(self, signum, frame):                                    
        raise TimeoutError(self.error_message)                                  
    def __enter__(self):                                                        
        signal.signal(signal.SIGALRM, self.handle_timeout)                      
        signal.alarm(self.seconds)                                              
    def __exit__(self, type, value, traceback):                                 
        signal.alarm(0) 

def check_targets():
    for folder in ['targets']:
        infolder = 'xyz_{}'.format(folder)
        outfile = 'res_{}.txt'.format(folder)
        xyzs = sorted(os.listdir(infolder))
        with open(outfile, 'w') as f:
            pass
        for infile in tqdm(xyzs):
            name = infile.split('.')[0]
            res = alexander('{}/{}'.format(infolder, infile), max_cross = 60, parallel_workers=False)
            try:
                is_knot = res['0_1'] < 0.5
            except KeyError:
                is_knot = True
            formatted_res = ','.join(['{}:{}'.format(k,v) for k,v in sorted(res.items())])
            is_knot_formatted = 'KNOT' if is_knot else '.   '
            #print(formatted_res)
            with open(outfile, 'a+') as f:
                f.write('{:<60} {} {}\n'.format(name, is_knot_formatted, formatted_res))

def check_models():
    for folder in ['models']:
        infolder = 'xyz_{}'.format(folder)
        outfile = 'res_{}.txt'.format(folder)
        xyzs = sorted(os.listdir(infolder))
        checked = set([])
        with open(outfile, 'r') as f:
            for line in f.readlines():
                model, num, *_ = line.split()
                checked.add('{}_{}'.format(model, num))
        # tqdm is not a necessary library
        for infile in tqdm(xyzs):
            name = infile.split('.')[0]
            if name in checked:
                continue
            name, num = name.split('_')
            try:
                with timeout(seconds=6):
                    # defaultly closure=2 which means two-point probabilistic closure,
                    # change to closure=0 for deterministic direct closure of endpoints.
                    # for some models max_cross needs to be increased
                    res = alexander('{}/{}'.format(infolder, infile), max_cross=200,
                                    closure=2, parallel_workers=False)
                try:
                    is_knot = res['0_1'] < 0.5
                except KeyError:
                    is_knot = True
                formatted_res = ','.join(['{}:{}'.format(k,v) for k,v in sorted(res.items())])
                is_knot_formatted = 'KNOT' if is_knot else '.   '
            except IndexError:
                is_knot_formatted = 'IE'
                formatted_res = 'IE'
            except TimeoutError:
                is_knot_formatted = 'TO'
                formatted_res = 'TO'
            except OverflowError:
                is_knot_formatted = 'OO'
                formatted_res = 'OO'
            with open(outfile, 'a+') as f:
                f.write('{:<14} {:2} {} {}\n'.format(name, num, is_knot_formatted, formatted_res))

if __name__ == '__main__':
    check_targets()
    check_models()
