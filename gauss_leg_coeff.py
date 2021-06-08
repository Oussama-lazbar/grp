from sys import argv
try :
    from bs4 import BeautifulSoup
    import requests as rq
except :
    print('Oups something went wrong. Check you have the packages needed : beautifulsoup4 and requests')
    exit
    
def GaussLegendre(n):
    """(wi,xi) de Gauss Legendre à n points. """
    #url = rq.get("https://pomax.github.io/bezierinfo/legendre-gauss.html")
    url = rq.get("https://web.archive.org/web/20210507045521/https://pomax.github.io/bezierinfo/legendre-gauss.html#n48")
    html = url.content

    soup = BeautifulSoup(html, 'html.parser')
    
    try :
        table = soup.find(id = 'n'+str(n)).find('tbody').find_all('td')
        with open(f'n{n}.txt', 'w') as f:
            for i in range(int(len(table)/3)):
                for j in range(1,3):
                    f.write(table[3*i  + j].text.strip())
                    f.write('   ')
                f.write('\n')
    except :
        print('This number of point is unreachable. Try with 1 < n < 65, or consider looking in another website.')
        exit


if __name__ == '__main__':
    GaussLegendre(int(argv[1]))
