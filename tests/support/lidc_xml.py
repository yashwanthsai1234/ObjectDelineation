"""Crafted LIDC read-message XML (objdel.data.lidc): ROIs, nodules, sessions, files."""


def square(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def roi_xml(z, points, inclusion=True):
    edges = "".join(f"<edgeMap><xCoord>{x}</xCoord><yCoord>{y}</yCoord></edgeMap>" for x, y in points)
    return (f"<roi><imageZposition>{z}</imageZposition><imageSOP_UID>s{z}</imageSOP_UID>"
            f"<inclusion>{'TRUE' if inclusion else 'FALSE'}</inclusion>{edges}</roi>")


def nodule_xml(rois, characteristics=True, nid="N"):
    chars = "<characteristics><subtlety>3</subtlety></characteristics>" if characteristics else ""
    rois = "".join(roi_xml(*r) for r in rois)
    return f"<unblindedReadNodule><noduleID>{nid}</noduleID>{chars}{rois}</unblindedReadNodule>"


def lidc_xml(sessions, series="1.2.3", root="LidcReadMessage"):
    body = "".join(f"<readingSession><servicingRadiologistID>anon</servicingRadiologistID>{''.join(s)}</readingSession>"
                   for s in sessions)
    return (f'<?xml version="1.0" encoding="UTF-8"?><{root} xmlns="http://www.nih.gov"><ResponseHeader>'
            f"<SeriesInstanceUid>{series}</SeriesInstanceUid></ResponseHeader>{body}</{root}>")


def write(tmp_path, sessions, name="001.xml", **kw):
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(lidc_xml(sessions, **kw))
    return p
