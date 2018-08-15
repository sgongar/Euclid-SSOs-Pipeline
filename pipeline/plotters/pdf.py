#!/usr/bin/python
# -*- coding: utf-8 -*-

from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas


def plot_pdf(logger, prefs_dict, result, full_cats, analysis_dict, conf_num):
    """

    @param logger:
    @param prefs_dict:
    @param result:
    @param full_cats:
    @param analysis_dict:
    @param conf_num:

    @return True: if everything goes alright
    """

    empty_SS0 = '/pcdisk/holly/sgongora/Documents/Euclid/empty_sso.png'
    logger.info('creating pdf file')
    x_page = [10, 10, 10, 10]
    y_page = [625, 420, 215, 10]

    c = canvas.Canvas(prefs_dict['report_out'] + '{}.pdf'.format(conf_num),
                      pagesize=A4)
    page_header_1 = 'detect_minarea    detect_thresh'
    page_header_2 = '    analysis_thresh          deblend_nthresh'
    page_header_3 = '    deblend_mincount'
    page_header = page_header_1 + page_header_2 + page_header_3
    c.setFont('Helvetica', 12)
    c.drawString(625, 10, page_header)

    p_1 = '{}'.format(analysis_dict['detect_minarea'])
    p_2 = '{}'.format(analysis_dict['detect_thresh'])
    p_3 = '{}'.format(analysis_dict['analysis_thresh'])
    p_4 = '{}'.format(analysis_dict['deblend_nthresh'])
    p_5 = '{}'.format(analysis_dict['deblend_mincount'])

    c.setFont('Helvetica', 10)
    for j in range(1, 6, 1):
        c.drawString(10 + (50 * j), 500, p_1)
        c.drawString(10 + (50 * j), 500, p_2)
        c.drawString(10 + (50 * j), 500, p_3)
        c.drawString(10 + (50 * j), 500, p_4)
        c.drawString(10 + (50 * j), 500, p_5)

    c.showPage()

    for source in range(0, len(result.index), 4):
        # four images per page
        for i in range(0, 4, 1):
            try:
                source_cat = int(result['source_cat'].iloc[source + i])
                source_num = int(result['source_num'].iloc[source + i])

                output_image = result['where_out'].iloc[source + i]
                logger.debug('placing image {}'.format(output_image))
                try:
                    c.drawImage(output_image, x_page[i], y_page[i],
                                width=288, height=175)
                except IOError:
                    logger.error('image {} not found'.format(output_image))
                    c.drawImage(empty_SS0, x_page[i], y_page[i],
                                width=288, height=175)

                x, y = x_page[i] + 320, y_page[i]

                logger.debug('writing legend for {}'.format(output_image))
                # write file name
                output_image = output_image.replace(prefs_dict['images_out'],
                                                    '')
                file_name = 'File: {}'.format(output_image)
                c.setFont('Helvetica', 10)
                c.drawString(x, y + 160, file_name)

                # write detections number
                detections_num = result['detections'].iloc[source + i]
                if detections_num != 1:
                    text_font = 'Helvetica-Bold'
                else:
                    text_font = 'Helvetica'
                detections_text = 'Detections: {}'.format(detections_num)

                c.setFont(text_font, 10)
                c.drawString(x, y + 140, detections_text)

                # write original values
                original_mag = int(result['original_mag'].iloc[source + i])
                original_pm = float(result['original_pm'].iloc[source + i])
                text_values = 'Magnitude: {} PM: {}'.format(original_mag,
                                                            original_pm)
                c.setFont('Helvetica', 10)
                c.drawString(x, y + 120, text_values)

                # write header
                table_header_1 = 'iso_area       threshold'
                table_header_2 = '    flags          mag_auto'
                table_header = table_header_1 + table_header_2
                c.setFont('Helvetica', 10)
                c.drawString(x, y + 100, table_header)

                # write data
                df_cat = full_cats[full_cats['source_cat'].isin([source_cat])]
                df_source = df_cat[df_cat['source_num'].isin([source_num])]

                iso_area_column = df_source['iso_area']
                threshold_column = df_source['threshold']
                flags_column = df_source['flags']
                mags_column = df_source['mag_auto']

                x_1 = x
                y_1 = y + 100
                for i in range(iso_area_column.size):
                    iso_area = iso_area_column.iloc[i]
                    threshold = '{0:.4f}'.format(threshold_column.iloc[i])
                    flags = flags_column.iloc[i]
                    mag_auto = '{0:.4f}'.format(mags_column.iloc[i])

                    t_1 = '{}'.format(iso_area)
                    t_2 = '{}'.format(threshold)
                    t_3 = '{}'.format(flags)
                    t_4 = '{}'.format(mag_auto)

                    c.setFont('Helvetica', 10)
                    c.drawString(x_1, y_1 - 20 * (i + 1), t_1)
                    c.drawString(x_1 + 60, y_1 - 20 * (i + 1), t_2)
                    c.drawString(x_1 + 110, y_1 - 20 * (i + 1), t_3)
                    c.drawString(x_1 + 160, y_1 - 20 * (i + 1), t_4)

            except IndexError:
                logger.debug('pdf finished')
        c.showPage()

    c.save()

    return True
